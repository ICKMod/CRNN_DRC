using CSV, DataFrames, PrettyTables, Printf, DifferentialEquations, DiffEqSensitivity, ForwardDiff, Plots, Markdown, NLsolve, SteadyStateDiffEq, LinearAlgebra
pyplot()
# --- Step 0: Pretty-print CSV (normalize reversible arrow to ⇌)

path = length(ARGS) > 0 ? ARGS[1] : "data.csv"
df = CSV.read(path, DataFrame; normalizenames=true)

const ANY_ARROW = r"(?:⇌|↔|<->|<=>|⇄|⟷|->|→|←)"

function _find_reaction_column(df)
    cols = collect(names(df))
    best_idx, best_score = 0, -1
    for (i, c) in pairs(cols)
        key  = c isa Symbol ? c : Symbol(c)
        vals = df[!, key]
        svals = (x -> (x isa AbstractString ? String(x) : string(x))).(vals)
        score = count(s -> occursin(ANY_ARROW, s), svals)
        if score > best_score
            best_idx, best_score = i, score
        end
    end
    best_idx == 0 && error("No reaction column detected.")
    return cols[best_idx]
end

rx_col_name = (:Elementary_Step in names(df)) ? :Elementary_Step :
              (:Step in names(df) ? :Step : _find_reaction_column(df))
rx_key = rx_col_name isa Symbol ? rx_col_name : Symbol(rx_col_name)

normalize_rev(s::AbstractString) = begin
    t = String(s)
    t = replace(t, r"\s*<\s*-\s*>\s*" => " ⇌ ")
    t = replace(t, r"\s*<\s*=\s*>\s*" => " ⇌ ")
    t = replace(t, "↔" => " ⇌ ")
    t = replace(t, "⟷" => " ⇌ ")
    t = replace(t, "⇄" => " ⇌ ")
    return t
end

df_disp = copy(df)
col = df[!, rx_key]
if any(ismissing, col)
    newcol = Vector{Union{Missing,String}}(undef, length(col))
    for i in eachindex(col)
        x = col[i]
        newcol[i] = ismissing(x) ? missing : normalize_rev(String(x))
    end
    df_disp[!, rx_key] = newcol
else
    df_disp[!, rx_key] = map(x -> normalize_rev(x), col)
end

num_idx = [i for (i, c) in enumerate(names(df_disp)) if eltype(df_disp[!, c]) <: Union{Missing, Real}]
fmt = (v, _i, j) -> (j in num_idx && v isa Real ? @sprintf("%.3e", v) : v)

pretty_table(df_disp;
    tf = tf_unicode_rounded,
    show_subheader = false,
    title = "DRC Input ($(basename(path)))",
    alignment = :l,
    show_row_number = true,
    row_number_column_title = "Step ID",
    formatters = fmt,
    crop = :none,
)

# ---------- Step 1: species extraction & classification ----------

rev_arrows = r"⇌|<->|↔|<=>|⇄|⟷"
fwd_arrows = r"->|→|←"

function find_reaction_column(df)
    cols = collect(names(df))                    
    best_idx, best_score = 0, -1
    for (i, c) in pairs(cols)
        key  = c isa Symbol ? c : Symbol(c)
        vals = df[!, key]
        svals = (x -> (x isa AbstractString ? String(x) : "")).(vals)
        score = count(s -> occursin(ANY_ARROW, s), svals)
        if score > best_score
            best_idx, best_score = i, score
        end
    end
    if best_idx == 0 || best_score == 0
        error("Could not auto-detect reaction column (no arrows found). Headers: $(names(df))")
    end
    return cols[best_idx]
end

step_col_name = (:Step in names(df)) ? :Step : find_reaction_column(df)
step_key = step_col_name isa Symbol ? step_col_name : Symbol(step_col_name)
println(">> Using reaction column: ", step_col_name)

function extract_species_from_step(step::AbstractString)::Vector{String}
    sides = split(step, ANY_ARROW)
    parts = Iterators.flatten(split.(strip.(sides), '+'))
    out = String[]
    for raw in parts
        tok = strip(raw)
        isempty(tok) && continue
        m  = match(r"^\s*(\d+(?:\.\d+)?)?\s*([^\s\+]+)\s*$", tok)
        sp = m === nothing ? tok : String(m.captures[2])

        occursin(r"^\d+(?:\.\d+)?(?:[eE][\+\-]?\d+)?$", sp) && continue
        occursin(r"[A-Za-z\*]", sp) || continue

        push!(out, sp)
    end
    return out
end

phase_of(sp::AbstractString) =
    sp == "*"         ? "free_site" :
    endswith(sp, "*") ? "adsorbed"  :
                        "gas"

raw_steps = df[!, step_key]
steps = String[]
for x in raw_steps
    ismissing(x) && continue
    s = x isa AbstractString ? String(x) : string(x)
    s = strip(s)
    isempty(s) && continue
    occursin(ANY_ARROW, s) || continue
    push!(steps, s)
end

is_reversible_base = [occursin(rev_arrows, s) for s in steps]

token_lists = [extract_species_from_step(s) for s in steps]
all_species = sort!(unique!(collect(Iterators.flatten(token_lists))))

catalog = DataFrame(Species = all_species, Phase = phase_of.(all_species))
phase_order = Dict("free_site"=>1, "adsorbed"=>2, "gas"=>3)
phase_rank(sp::AbstractString) = get(phase_order, sp, 99)
catalog[!, :_rank] = phase_rank.(catalog.Phase)

sort!(catalog, [:_rank, :Species])
select!(catalog, [:Species, :Phase])


try
    pretty_table(catalog;
        header = ["Species", "Phase"],
        show_subheader = false,
        tf = tf_unicode_rounded,
        title = "Species Catalog (by phase) — column: $(string(step_col_name))",
        alignment = :l,
        show_row_number = true,
        row_number_column_title = "#",
        crop = :none,
    )
catch
    pretty_table(catalog;
        header = ["Species", "Phase"],
        show_subheader = false,
        tf = tf_unicode_rounded,
        title = "Species Catalog (by phase) — column: $(string(step_col_name))",
        crop = :none,
    )
end

# --- Step 2: Stoichiometric matrix for adsorbates + free site (mixed rev/irrev)

split_arrow(s::AbstractString) = begin
    parts = split(s, ANY_ARROW)
    length(parts) >= 2 ? (strip(parts[1]), strip(parts[end])) : (strip(s), "")
end

function parse_side(side::AbstractString)
    d = Dict{String,Float64}()
    for tok in split(side, '+')
        t = strip(tok); isempty(t) && continue
        m = match(r"^\s*(\d+(?:\.\d+)?)?\s*([^\s\+]+)\s*$", t)
        if m === nothing
            sp, ν = t, 1.0
        else
            ν = m.captures[1] === nothing ? 1.0 : parse(Float64, m.captures[1])
            sp = String(m.captures[2])
        end
        occursin(r"^\d+(?:\.\d+)?(?:[eE][\+\-]?\d+)?$", sp) && continue
        d[sp] = get(d, sp, 0.0) + ν
    end
    return d
end

reactions_base = Vector{Tuple{Dict{String,Float64}, Dict{String,Float64}}}()
for s in steps
    lhs, rhs = split_arrow(s)
    L = parse_side(lhs); R = parse_side(rhs)
    push!(reactions_base, (L, R))
end
isempty(reactions_base) && error("Parsed 0 reactions from 'steps' — check arrow detection or input.")

reactions = Vector{Tuple{Dict{String,Float64}, Dict{String,Float64}}}()
base_id  = Int[]
dir_flag = String[]

for (j, (L, R)) in enumerate(reactions_base)
    push!(reactions, (L, R))
    push!(base_id, j)
    push!(dir_flag, "fwd")

    if is_reversible_base[j]
        push!(reactions, (R, L))
        push!(base_id, j)
        push!(dir_flag, "rev")
    end
end

non_gas_df = filter(:Phase => p -> p != "gas", catalog)
cols_species = String.(non_gas_df.Species)

if any(haskey(L, "*") || haskey(R, "*") for (L, R) in reactions) && !("*" in cols_species)
    cols_species = vcat("*", cols_species)
end
if "*" in cols_species
    cols_species = vcat("*", filter(!=("*"), cols_species))
end

S_ads = zeros(Float64, length(reactions), length(cols_species))
for (j, (L, R)) in enumerate(reactions)
    for (i, sp) in pairs(cols_species)
        S_ads[j, i] = get(R, sp, 0.0) - get(L, sp, 0.0)
    end
end

display_species = String.(catalog.Species)

if any(haskey(L, "*") || haskey(R, "*") for (L, R) in reactions) && !("*" in display_species)
    display_species = vcat("*", display_species)
end
if "*" in display_species
    display_species = vcat("*", filter(!=("*"), display_species))
end

S_disp = zeros(Float64, length(reactions), length(display_species))
for (j, (L, R)) in enumerate(reactions)
    for (i, sp) in pairs(display_species)
        S_disp[j, i] = get(R, sp, 0.0) - get(L, sp, 0.0)
    end
end

header_vec = vcat(["BaseRx"], display_species)
disp = Array{Any}(undef, size(S_disp, 1), 1 + length(display_species))

for j in eachindex(base_id)
    bid = base_id[j]
    is_rev = is_reversible_base[bid]
    if is_rev
        disp[j, 1] = string(bid, dir_flag[j] == "rev" ? "r" : "f")
    else
        disp[j, 1] = string(bid)
    end
end

for j in 1:length(display_species)
    disp[:, 1 + j] = S_disp[:, j]
end

pretty_table(disp;
    header = header_vec,
    show_subheader = false,
    tf = tf_unicode_rounded,
    title = "Stoichiometric matrix",
    alignment = :r,
    crop = :none,
)

# --- Step 2b: raw matrix + row-major flatten
S_ads_numeric = copy(S_ads)
S_ads_flat_rowmajor = vec(permutedims(S_ads_numeric))

# --- Step 3: Interactive inputs for non-gas coverages and gas pressures

non_gas_df = filter(:Phase => p -> p != "gas", catalog)
gas_df     = filter(:Phase => p -> p == "gas", catalog)
non_gas_species = collect(non_gas_df.Species)
gas_species     = collect(gas_df.Species)

function _eval_ast(ex)
    ex isa Number && return float(ex)
    ex isa Symbol && error("symbols not allowed")
    ex isa Expr || error("unsupported")
    ex.head == :call || error("unsupported")
    op = ex.args[1]; args = ex.args[2:end]; vals = map(_eval_ast, args)
    op === :+ && return +(vals...)
    op === :- && return (length(vals)==1 ? -vals[1] : -(vals...))
    op === :* && return *(vals...)
    op === :/ && return /(vals...)
    op === :^ && return ^(vals...)
    error("operator $op not allowed")
end
safe_eval(s::AbstractString) = _eval_ast(Meta.parse(s))


function prompt_values(names::Vector{String}, label::String)
    println("Enter $label for $(length(names)) species (number or arithmetic expression).")
    vals = Float64[]; raws = String[]
    for (i, sp) in enumerate(names)
        while true
            try
                print("[$i/$(length(names))] $sp = ")
                inp = strip(readline())
                isempty(inp) && (println("  empty input, try again."); continue)
                v = safe_eval(inp)
                push!(vals, v); push!(raws, inp)
                break
            catch e
                println("  invalid entry: $e")
            end
        end
    end
    return vals, raws
end

prs_vals, prs_raws = prompt_values(gas_species, "pressure")

gas_input = DataFrame(
    "Species" => gas_species,
    "Phase" => gas_df.Phase,
    "Expression" => prs_raws,
    "Pressure" => prs_vals
)

pretty_table(gas_input; header=String.(names(gas_input)), show_subheader=false, 
             tf=tf_unicode_rounded, title="Gas inputs (pressure)", alignment=:l, crop=:none)

pressure_map = Dict(gas_species[i] => prs_vals[i] for i in eachindex(gas_species))

println(pressure_map)

function ask_steady_state()
    ss_choice = lowercase(strip(_read_choice(
        "Go to steady state? (y/n)", 
        "n"; 
        envkey="DRC_STEADY_STATE")))
    return ss_choice in ["y", "yes", "true", "1"]
end

const R_gas = 8.314462618

const _LAST_T_ARR     = Ref{Union{Nothing,Float64}}(nothing)
const _LAST_EA_UNIT   = Ref{Union{Nothing,String}}(nothing)
const _LAST_B_ARR     = Ref{Union{Nothing,Float64}}(nothing)
const _LAST_EAPP_UNIT = Ref{Union{Nothing,String}}(nothing)

canon_case(s) = replace(String(s), r"[^A-Za-z0-9]" => "")

function find_col_case(df, candidates::Vector{String})
    targets = Set(canon_case(c) for c in candidates)
    for nm in names(df)
        if canon_case(nm) in targets
            return nm isa Symbol ? nm : Symbol(nm)
        end
    end
    error("Could not find any of: $(join(candidates, ", ")). Headers: $(names(df))")
end

function has_col_case(df, candidates::Vector{String})::Bool
    targets = Set(canon_case(c) for c in candidates)
    for nm in names(df)
        if canon_case(nm) in targets
            return true
        end
    end
    return false
end

_to_float(x) = x isa Real ? float(x) : parse(Float64, replace(String(x), 'd' => 'e'))

normalize_K(x) = begin
    if x === missing || x == ""
        return Inf
    elseif x isa AbstractString
        s = strip(lowercase(String(x)))
        s in ("inf", "infinity", "∞") && return Inf
        return parse(Float64, replace(s, 'd'=>'e'))
    else
        return Float64(x)
    end
end

kr_from_kfK(kf::Real, K::Real) = isfinite(K) ? (kf / K) : 0.0

_read_str(prompt, default) = begin
    print("$(prompt) [default: $(default)]: "); flush(stdout)
    local s::String
    try
        s = strip(readline(stdin))
        isempty(s) ? string(default) : s
    catch
        string(default)
    end
end

function _read_float(prompt::AbstractString, default::Real; envkey::Union{Nothing,String}=nothing)
    key = envkey === nothing ? replace(uppercase(prompt), " " => "_") : envkey
    val = get(ENV, key, "")
    if !isempty(val)
        try return parse(Float64, val) catch end
    end
    parse(Float64, _read_str(prompt, default))
end

_read_choice(prompt::AbstractString, default::AbstractString; envkey::Union{Nothing,String}=nothing) = begin
    key = envkey === nothing ? replace(uppercase(prompt), " " => "_") : envkey
    val = get(ENV, key, "")
    if !isempty(val)
        return String(val)
    end
    _read_str(prompt, default)
end

function arrhenius_k(A, Ea, T; Ea_unit::AbstractString="kJ/mol", b::Real=0.0)
    EaJ = lowercase(Ea_unit) in ("kj/mol","kjmol","kj") ? Ea*1000.0 : Ea
    return A * (T^b) * exp(-EaJ / (R_gas * T))
end

function _thermo_check_logtol()
    val = strip(get(ENV, "DRC_THERMO_CHECK_LOGTOL", "1e-6"))
    try
        return parse(Float64, val)
    catch
        @warn "Could not parse DRC_THERMO_CHECK_LOGTOL='$val'; using 1e-6."
        return 1e-6
    end
end

function _check_K_kf_kr_consistency(row_i::Int, reaction_string::AbstractString,
                                   K_given::Real, kf_r::Real, kr_r::Real;
                                   tolerance_logK::Float64=_thermo_check_logtol())
    if !isfinite(K_given)
        if kr_r == 0.0 || !isfinite(kf_r / kr_r)
            return true
        end
        error("Thermodynamic consistency check failed at row $row_i: provided K=$(K_given) but kf/kr=$(kf_r/kr_r). Check K, kf, and kr. Reaction: $(reaction_string)")
    end
    K_given <= 0.0 && error("Equilibrium constant K must be positive at row $row_i. Reaction: $(reaction_string)")
    K_rates = kr_r == 0.0 ? Inf : (kf_r / kr_r)
    if !(isfinite(K_rates) && K_rates > 0.0)
        error("Thermodynamic consistency check failed at row $row_i: provided K=$(K_given) but kf/kr=$(K_rates). Check K, kf, and kr. Reaction: $(reaction_string)")
    end
    log_mismatch = abs(log(K_rates) - log(K_given))
    if log_mismatch > tolerance_logK
        error("Thermodynamic consistency check failed at row $row_i: provided K=$(K_given) but kf/kr=$(K_rates). Check K, kf, and kr. kf=$(kf_r), kr=$(kr_r), log mismatch=$(log_mismatch). Reaction: $(reaction_string)")
    end
    return true
end

function _get_arr_cols(df; which::Symbol)
    if which === :fwd
        Akey  = find_col_case(df, ["Af","A_f","Aforward","AForward","A_fwd","A_F"])
        EAkey = find_col_case(df, ["Eaf","Ea_f","EaForward","Ea_forward","Eafwd","Ea_F"])
    elseif which === :rev
        Akey  = find_col_case(df, ["Ar","A_r","Areverse","AReverse","A_rev","A_R"])
        EAkey = find_col_case(df, ["Ear","Ea_r","EaReverse","Ea_reverse","Earv","Ea_R"])
    else
        error("which must be :fwd or :rev")
    end
    return Akey, EAkey
end

# --- Step 4: Allocate ODE vector (using non-gas species only)

use_steady_state = ask_steady_state()

ss_converged = false
t_ss = 0.0
ss_sol = nothing
u_ss = nothing

function build_rate_constants(df, rx_key, reactions, is_reversible_base)
    println("=== Kinetics data mode (auto-detected) ===")

    pair_env = lowercase(strip(get(ENV, "DRC_PAIR", "")))
    src_env  = lowercase(strip(get(ENV, "DRC_SOURCE", "")))
    has_K  = has_col_case(df, ["K","Keq","K_eq","KEquilibrium","Kconst","Keqconst"])
    has_kr = has_col_case(df, ["kr","k_r","kreverse","kReverse"])
    has_rev_arr = has_col_case(df, ["Ar","A_r","Areverse","AReverse","A_rev","A_R"]) ||
                  has_col_case(df, ["Ear","Ea_r","EaReverse","Ea_reverse","Earv","Ea_R"])

    src_choice = ""
    if !isempty(src_env)
        src_choice = src_env
    else
        has_Af  = has_col_case(df, ["Af","A_f","Aforward","AForward","A_fwd","A_F"])
        has_Eaf = has_col_case(df, ["Eaf","Ea_f","EaForward","Ea_forward","Eafwd","Ea_F"])
        has_Ar  = has_col_case(df, ["Ar","A_r","Areverse","AReverse","A_rev","A_R"])
        has_Ear = has_col_case(df, ["Ear","Ea_r","EaReverse","Ea_reverse","Earv","Ea_R"])
        src_choice = (has_Af || has_Eaf || has_Ar || has_Ear) ? "arrhenius" : "direct"
    end
    src_choice ∈ ["direct","arrhenius"] || error("Parameter source must be 'Direct' or 'Arrhenius'")
    thermo_check_required = has_K && (
        (src_choice == "direct" && has_kr) ||
        (src_choice == "arrhenius" && has_rev_arr)
    )
    reverse_input_available = src_choice == "direct" ? has_kr : has_rev_arr

    pair_choice = ""
    if !isempty(pair_env)
        pair_choice = pair_env
    else
        if thermo_check_required
            pair_choice = "(keq,kf)"
        elseif has_K && !reverse_input_available
            pair_choice = "(keq,kf)"
        elseif !has_K && reverse_input_available
            pair_choice = "(kf,kr)"
        else
            error("Could not infer kinetic pair type from headers. Expected either (K/Keq with k/kf) or (kf with kr) or Arrhenius (Af/Eaf with optional Ar/Ear). Headers: $(names(df))")
        end
    end

    pair_choice ∈ ["(keq,kf)", "(kf,kr)"] || error("Kinetic Constants pair type must be '(Keq,kf)' or '(kf,kr)'")

    println(">> Pair type: ", pair_choice == "(keq,kf)" ? "(Keq,kf)" : "(kf,kr)")
    println(">> Parameter source: ", src_choice == "direct" ? "Direct" : "Arrhenius")

    valid_rows = Int[]
    for r in 1:nrow(df)
        s_any = df[r, rx_key]
        s = s_any isa AbstractString ? String(s_any) : string(s_any)
        occursin(ANY_ARROW, s) && push!(valid_rows, r)
    end
    length(valid_rows) == length(reactions) || @warn "Row count for reactions differs from parsed reactions; proceeding by order."

    kf = Float64[]; Keq = Float64[]
    checked_rows = 0

    if src_choice == "direct"
        if pair_choice == "(keq,kf)"
            k_key = find_col_case(df, ["k","kf","k_f","kforward","kForward"])
            K_key = find_col_case(df, ["K","Keq","K_eq","KEquilibrium","Kconst","Keqconst"])
            kr_key = thermo_check_required ? find_col_case(df, ["kr","k_r","kreverse","kReverse"]) : nothing
            for (j, r) in enumerate(valid_rows)
                kf_r = _to_float(df[r, k_key])
                K_given = normalize_K(df[r, K_key])
                if thermo_check_required && is_reversible_base[j]
                    kr_r = _to_float(df[r, kr_key])
                    _check_K_kf_kr_consistency(r, string(df[r, rx_key]), K_given, kf_r, kr_r)
                    checked_rows += 1
                end
                push!(kf, kf_r)
                push!(Keq, K_given)
                @assert Keq[end] != 0.0 "Equilibrium constant K must not be 0 at row $r."
            end
        else
            kf_key = find_col_case(df, ["kf","k","k_f","kforward","kForward"])
            kr_key = find_col_case(df, ["kr","k_r","kreverse","kReverse"])
            K_key = thermo_check_required ? find_col_case(df, ["K","Keq","K_eq","KEquilibrium","Kconst","Keqconst"]) : nothing
            for (j, r) in enumerate(valid_rows)
                kf_r = _to_float(df[r, kf_key])
                kr_r = _to_float(df[r, kr_key])
                if thermo_check_required && is_reversible_base[j]
                    K_given = normalize_K(df[r, K_key])
                    _check_K_kf_kr_consistency(r, string(df[r, rx_key]), K_given, kf_r, kr_r)
                    checked_rows += 1
                end
                push!(kf, kf_r)
                push!(Keq, kr_r == 0.0 ? Inf : (kf_r/kr_r))
            end
        end
    else
        if _LAST_T_ARR[] === nothing
            _LAST_T_ARR[] = _read_float("Temperature (K)", 500.0)
        end
        if _LAST_EA_UNIT[] === nothing
            _LAST_EA_UNIT[] = _read_choice("Ea unit (kJ/mol or J/mol)", "kJ/mol"; envkey="DRC_EA_UNIT")
        end
        if _LAST_B_ARR[] === nothing
            _LAST_B_ARR[] = _read_float("Arrhenius exponent b", 0.0; envkey="DRC_B_ARR")
        end
        T = _LAST_T_ARR[]
        Ea_unit = _LAST_EA_UNIT[]
        b = _LAST_B_ARR[]

        if pair_choice == "(keq,kf)"
            Af_key, Eaf_key = _get_arr_cols(df; which=:fwd)
            K_key = find_col_case(df, ["K","Keq","K_eq","KEquilibrium","Kconst","Keqconst"])
            rev_keys = thermo_check_required ? _get_arr_cols(df; which=:rev) : (nothing, nothing)
            Ar_key, Ear_key = rev_keys
            for (j, r) in enumerate(valid_rows)
                Af  = _to_float(df[r, Af_key])
                Eaf = _to_float(df[r, Eaf_key])
                kf_r = arrhenius_k(Af, Eaf, T; Ea_unit=Ea_unit, b=b)
                K_given = normalize_K(df[r, K_key])
                if thermo_check_required && is_reversible_base[j]
                    Ar  = _to_float(df[r, Ar_key])
                    Ear = _to_float(df[r, Ear_key])
                    kr_r = arrhenius_k(Ar, Ear, T; Ea_unit=Ea_unit, b=b)
                    _check_K_kf_kr_consistency(r, string(df[r, rx_key]), K_given, kf_r, kr_r)
                    checked_rows += 1
                end
                push!(kf, kf_r)
                push!(Keq, K_given)
                @assert Keq[end] != 0.0 "Equilibrium constant K must not be 0 at row $r."
            end
        else
            Af_key, Eaf_key = _get_arr_cols(df; which=:fwd)
            Ar_key, Ear_key = _get_arr_cols(df; which=:rev)
            K_key = thermo_check_required ? find_col_case(df, ["K","Keq","K_eq","KEquilibrium","Kconst","Keqconst"]) : nothing
            for (j, r) in enumerate(valid_rows)
                Af  = _to_float(df[r, Af_key])
                Eaf = _to_float(df[r, Eaf_key])
                Ar  = _to_float(df[r, Ar_key])
                Ear = _to_float(df[r, Ear_key])
                kf_r = arrhenius_k(Af, Eaf, T; Ea_unit=Ea_unit, b=b)
                kr_r = arrhenius_k(Ar, Ear, T; Ea_unit=Ea_unit, b=b)
                if thermo_check_required && is_reversible_base[j]
                    K_given = normalize_K(df[r, K_key])
                    _check_K_kf_kr_consistency(r, string(df[r, rx_key]), K_given, kf_r, kr_r)
                    checked_rows += 1
                end
                push!(kf, kf_r)
                push!(Keq, kr_r == 0.0 ? Inf : (kf_r/kr_r))
            end
        end
    end

    if thermo_check_required && checked_rows > 0
        println(">> Thermodynamic consistency check passed for rows with K, kf, and kr.")
    end

    @assert length(kf)  == length(reactions) "kf rows ≠ reaction count."
    @assert length(Keq) == length(reactions) "K rows ≠ reaction count."

    return kf, Keq
end

@assert isdefined(Main, :reactions) "Step 2 reactions not found."

k_fwd_raw, K_eq_raw = build_rate_constants(df, rx_key, reactions_base, is_reversible_base)

for j in eachindex(K_eq_raw)
    if !is_reversible_base[j]
        K_eq_raw[j] = Inf
    end
end

k_rev_raw = map(kr_from_kfK, k_fwd_raw, K_eq_raw)
k_fwd_eff = copy(k_fwd_raw)
k_rev_eff = copy(k_rev_raw)

@assert isdefined(Main, :catalog) "Step 1 catalog not found."
@assert isdefined(Main, :pressure_map) "Step 3 pressure_map not found."

gas_set = Set(String.(filter(:Phase => ==("gas"), catalog).Species))

for (j, (L, R)) in enumerate(reactions_base)
    pf = 1.0
    for (sp, ν) in L
        if sp in gas_set
            @assert haskey(pressure_map, sp) "Missing pressure for gas species $sp."
            pf *= pressure_map[sp]^ν
        end
    end
    k_fwd_eff[j] *= pf

    pr = 1.0
    for (sp, ν) in R
        if sp in gas_set
            @assert haskey(pressure_map, sp) "Missing pressure for gas species $sp."
            pr *= pressure_map[sp]^ν
        end
    end
    k_rev_eff[j] *= pr
end

ode_rates = Float64[]
for j in eachindex(reactions_base)
    push!(ode_rates, k_fwd_eff[j])
    if is_reversible_base[j]
        push!(ode_rates, k_rev_eff[j])
    end
end

@assert length(ode_rates) == length(reactions)

n_rx_total = length(reactions)
n_non_gas  = nrow(non_gas_df)
ode_vec    = zeros(Float64, n_rx_total * (n_non_gas + 1))

@assert length(ode_vec) >= n_rx_total "ode_vec too small for rate constants."
ode_vec[1:n_rx_total] = log.(ode_rates)

@assert isdefined(Main, :S_ads_flat_rowmajor) "S_ads_flat_rowmajor not found (from Step 2b)."
n_cols = length(S_ads_flat_rowmajor) / n_rx_total
@assert n_rx_total * (n_cols + 1) == length(ode_vec) "ode_vec size mismatch."

start = n_rx_total + 1
ode_vec[start : start + length(S_ads_flat_rowmajor) - 1] = S_ads_flat_rowmajor;

function _read_int(prompt::AbstractString, default::Integer)
    val = get(ENV, replace(uppercase(prompt), " " => "_"), "")
    if !isempty(val)
        try return parse(Int, val) catch end
    end
    parse(Int, _read_str(prompt, default))
end

function _choose_solver(name::AbstractString)
    n = lowercase(strip(String(name)))
    table = Dict(
        "kvaerno5" => Kvaerno5(),
        "tsit5"    => Tsit5(),
        "vern9"    => Vern9(),
        "rk4"      => RK4(),
        "rodas5"   => Rodas5(),
        "rodas5p"  => Rodas5P(),
        "trbdf2"   => TRBDF2(),
        "rodas4p"  => Rodas4P()
    )
    get(table, n, Rodas5P())
end

function _default_ode_maxiters()
    env_val = strip(get(ENV, "DRC_MAXITERS", ""))
    if !isempty(env_val)
        parsed = try parse(Int, env_val) catch; nothing end
        if parsed !== nothing && parsed > 0
            return parsed
        end
    end
    return 200000
end

function get_ode_user_settings()
    println("\n=== ODE Settings ===")
    println("Solver reference: https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/\n")

    t0 = _read_float("Time Span: t0", 0.0)
    t1 = _read_float("Time Span: t1", 50.0)
    t1 <= t0 && error("t1 must be greater than t0")
    tspan = (t0, t1)

    println("\n--- Time grid for saving the profiles ---")
    save0  = _read_float("save_start", t0)
    save1  = _read_float("save_end",   t1)
    n_main = _read_int("save_points",  500)
    n_main < 2 && error("save_points must be ≥ 2")
    tsteps_main = range(save0, save1; length=n_main)

    default_alg = "Rodas5P"
    alg_name = get(ENV, "DRC_ALG", _read_str("solver (e.g., Tsit5, Vern9, Rodas5, Kvaerno5)", default_alg))
    alg = _choose_solver(alg_name)
    maxiters = _read_int("maxiters", _default_ode_maxiters())

    println("\nChosen settings:")
    println("  tspan       = $(tspan)")
    println("  save grid   = $(first(tsteps_main)) → $(last(tsteps_main))  (n=$(length(tsteps_main)))")
    println("  solver      = $(typeof(alg))")
    println("  maxiters    = $(maxiters)")

    return tspan, collect(tsteps_main), alg, maxiters
end

# --- Step 5: Define & solve ODEs, save solution

p = ode_vec

nr = length(reactions)
ns = size(S_ads_numeric, 2)
st_input = S_ads_numeric
lb, ub = 1e-100, 1e0

function p2vec(p)
    nr_loc = length(reactions)
    ns_loc = Int((length(p) - nr_loc) ÷ nr_loc)
    @assert nr_loc * (ns_loc + 1) == length(p) "Parameter vector length mismatch with reactions/stoichiometry"
    w_b   = @view p[1:nr_loc]
    w_out = reshape(@view(p[nr_loc+1:end]), ns_loc, nr_loc)
    w_in  = clamp.(-w_out, 0, 5)
    return w_in, w_b, w_out
end

function display_p(p)
    w_in, w_b, w_out = p2vec(p);
    println("species (column) reaction (row)")
    println("w_in")
    show(stdout, "text/plain", round.(w_in', digits=3))

    println("\nw_b")
    show(stdout, "text/plain", round.(exp.(w_b'), digits=6))

    println("\nw_out")
    show(stdout, "text/plain", round.(w_out', digits=3))
    println("\n\n")
end

function compute_reaction_rates(u, p)
    w_in, w_b, w_out = p2vec(p)
    @assert length(u) == size(w_out, 1) "State vector length $(length(u)) does not match stoichiometric columns $(size(w_out,1))"
    w_in_x = w_in' * @. log(clamp(u, lb, ub))
    rate_all = @. exp(w_in_x + w_b)
    return rate_all, w_out
end

function crnn!(du, u, p, t)
    rate_all, w_out = compute_reaction_rates(u, p)
    du .= w_out * rate_all;
end

@inline function _sigmoid(x)
    x > 40  && return 1.0
    x < -40 && return 0.0
    return 1.0 / (1.0 + exp(-x))
end

@inline function _logit(y)
    return log(y / (1.0 - y))
end

function _warm_start_guess(u0_init, p, lb, ub;
                           t_end = 1e3,
                           abstol = 1e-12,
                           reltol = 1e-10,
                           maxiters = _default_ode_maxiters())
    u0c = clamp.(u0_init, lb, ub)
    prob_ss = ODEProblem(crnn!, u0c, (0.0, t_end), p)

    sol = solve(prob_ss, Rodas5P();
                abstol=abstol, reltol=reltol,
                maxiters=maxiters,
                saveat=[t_end],
                isoutofdomain = (u, p, t) -> any(u .< lb))

    uw = Array(sol)[:, end]
    return clamp.(uw, lb, ub)
end

function solve_steady_state(u0_init, p, lb, ub;
                            warm_start=false,
                            t_end = 1e3,
                            ftol=1e-12,
                            maxiters=10000,
                            accept_resid=1e-6)

    function residual!(F, u)
        rate_all, w_out = compute_reaction_rates(u, p)
        F .= w_out * rate_all
    end

    u_guess = warm_start ? _warm_start_guess(u0_init, p, lb, ub; t_end=t_end) : clamp.(u0_init, lb, ub)

    F0 = similar(u_guess)
    residual!(F0, u_guess)
    scale = max.(abs.(F0), 1.0)

    function u_from_z(z)
        u = similar(z)
        @inbounds for i in eachindex(z)
            u[i] = lb + (ub - lb) * _sigmoid(z[i])
        end
        return u
    end

    function z_from_u(u)
        z = similar(u)
        @inbounds for i in eachindex(u)
            y = (u[i] - lb) / (ub - lb)
            y = clamp(y, 1e-14, 1.0 - 1e-14)
            z[i] = _logit(y)
        end
        return z
    end

    z0 = z_from_u(u_guess)

    function residual_z!(Fz, z)
        u = u_from_z(z)
        residual!(Fz, u)
        Fz ./= scale
        return nothing
    end

    println("\n=== Solving for Steady State (warm_start=$(warm_start), t_end=$(t_end)) ===")
    println("Initial guess: ")

    u_guess_df = DataFrame(
    Index = collect(1:length(non_gas_species)),
    Species = non_gas_species,
    u_guess = u_guess
    )

    pretty_table(u_guess_df;
    header=["Index", "Species", "u0 initial guess"],
    show_subheader=false,
    tf=tf_unicode_rounded,
    crop=:none)


    methods_to_try = [
        (:trust_region, "Trust Region"),
        (:anderson,     "Anderson Acceleration"),
        (:newton,       "Newton"),
        (:broyden,      "Broyden's quasi-Newton")
    ]

    best_u = u_guess
    best_res = Inf

    for (method, name) in methods_to_try
        println("Attempting $name method...")
        try
            result = nlsolve(residual_z!, z0;
                             method=method,
                             ftol=ftol,
                             iterations=maxiters,
                             show_trace=false)

            ok = NLsolve.converged(result)
            u_ss = u_from_z(result.zero)

            F_test = similar(u_ss)
            residual!(F_test, u_ss)
            max_residual = maximum(abs.(F_test))

            if ok
                println("✓ Converged with $name method. Max residual = $(max_residual)")
            else
                println("✗ Did not converge with $name method. Max residual = $(max_residual)")
            end

            if max_residual < best_res
                best_res = max_residual
                best_u = u_ss
            end

            if max_residual < accept_resid
                println("✓ Accepted based on residual threshold (ok=$(ok)).")
                return u_ss, true
            end
        catch e
            println("✗ $name method failed with error: $e")
        end
    end

    println("⚠ Could not find steady state with the provided initial guess/bounds. Best max residual = $(best_res)")
    return best_u, false
end

function solve_steady_state_dynamicss(u0_init, p, lb, ub;
    t_end::Real = 1e3,
    accept_resid::Real = 1e-6,
    alg = Rodas5P(),
    odesolve_abstol::Real = 1e-12,
    odesolve_reltol::Real = 1e-10,
    dtmin::Real = 1e-25,
    maxiters::Int = _default_ode_maxiters())

    u0c = clamp.(u0_init, lb, ub)

    prob_ss = ODEProblem(crnn!, u0c, (0.0, float(t_end)), p)

    cb = TerminateSteadyState(float(accept_resid), 0.0; min_t=1e-12)

    sol = solve(prob_ss, alg;
        callback = cb,
        abstol   = float(odesolve_abstol),
        reltol   = float(odesolve_reltol),
        dtmin    = float(dtmin),
        maxiters = Int(maxiters),
        save_everystep = false)

    u_ss = clamp.(sol.u[end], lb, ub)
    t_ss = sol.t[end]

    du = similar(u_ss)
    crnn!(du, u_ss, p, t_ss)
    max_residual = maximum(abs.(du))

    converged = max_residual < accept_resid
    return u_ss, converged, t_ss, sol
end

if use_steady_state
    println("\n=== Steady State Mode ===")
    println("Attempting to find steady state coverages...")
    
    ode_maxiters = _read_int("steady_state_maxiters", _default_ode_maxiters())
    ss_profile_points = _read_int("steady_state_profile_points", 500)
    ss_profile_points < 2 && (ss_profile_points = 2)
    ss_accept_resid = _read_float("steady_state_accept_threshold", 1e-8)
    ss_accept_resid <= 0 && error("steady_state_accept_threshold must be > 0")

    u0_guess = zeros(Float64, length(non_gas_species))
    free_site_idx = findfirst(==("*"), String.(non_gas_species))
    if free_site_idx !== nothing
        u0_guess[free_site_idx] = 1.0
    else
        @warn "Free-site '*' not found in non_gas_species; using first species as 1.0 for u0_guess."
        u0_guess[1] = 1.0
    end

    #u0_guess = Float64[0.68, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
    
    u_ss, ss_converged, t_ss, ss_sol = solve_steady_state_dynamicss(u0_guess, p, lb, ub;
                                                    t_end=1e3,
                                                    accept_resid=ss_accept_resid,
                                                    maxiters=ode_maxiters)
    
    if ss_converged
        println("\n✓ Steady state found successfully!")
        println("Steady state reached at t ≈ ", t_ss, " s")
        
        coverage_map = Dict(non_gas_species[i] => u_ss[i] for i in eachindex(non_gas_species))
        
        ss_df = DataFrame(
            Species = non_gas_species,
            SteadyState_Coverage = u_ss
        )
        pretty_table(ss_df; 
            header=["Species", "Steady State θ"],
            show_subheader=false, 
            tf=tf_unicode_rounded, 
            title="Steady State Coverages (computed)",
            formatters = (v, i, j) -> (j > 1 && v isa Real ? @sprintf("%.6e", v) : v),
            alignment=:r,
            crop=:none)
        
        # Use the trajectory from the steady-state solve for profiles; initial = first point
        u0 = ss_sol.u[1]
        cov_vals = u0
        cov_raws = [@sprintf("%.6e", v) for v in u0]
        
    else
        println("\n⚠ Steady state not found. Switching to transient mode.")
        println("Please provide initial coverages manually.\n")
        
        cov_vals, cov_raws = prompt_values(non_gas_species, "initial coverage")
        u0 = cov_vals
        coverage_map = Dict(non_gas_species[i] => cov_vals[i] for i in eachindex(non_gas_species))
    end
    
else
    println("\n=== Transient Mode ===")
    cov_vals, cov_raws = prompt_values(non_gas_species, "initial coverage")
    u0 = cov_vals
    coverage_map = Dict(non_gas_species[i] => cov_vals[i] for i in eachindex(non_gas_species))
end

non_gas_input = DataFrame(
    "Species" => non_gas_species,
    "Phase" => non_gas_df.Phase,
    "Expression" => cov_raws,
    "Initial Coverage" => cov_vals
)

pretty_table(non_gas_input; header=String.(names(non_gas_input)), show_subheader=false, 
             tf=tf_unicode_rounded, title="Non-gas inputs (coverage)", alignment=:l, crop=:none)

println(coverage_map)

if use_steady_state && ss_converged
    tspan = (0.0, t_ss)
    tsteps = collect(range(0.0, t_ss; length=ss_profile_points))
    alg = Rodas5P()
else
    tspan, tsteps, alg, ode_maxiters = get_ode_user_settings()
end

prob = ODEProblem(crnn!, u0, tspan, p)

function predict_neuralode(prob, u0, tsteps; _alg=Rodas5P(), _maxiters=_default_ode_maxiters())
    sol = solve(prob, _alg;
        u0=u0,
        saveat=tsteps,
        sensealg=ForwardDiffSensitivity(),
        abstol=1e-12,
        reltol=1e-12,
        dtmin=1e-25,
        maxiters=Int(_maxiters))

    if sol.retcode == ReturnCode.MaxIters
        @warn "Solve reached maxiters. Increase maxiters (or use a different stiff solver) for a longer trajectory."
    elseif sol.retcode !== ReturnCode.Success && sol.retcode !== ReturnCode.Terminated
        @warn "ODE solve ended with retcode $(sol.retcode). The profile grid may be truncated."
    end

    return clamp.(Array(sol), lb, ub), sol.t
end

sol_arr, out_t = predict_neuralode(prob, u0, tsteps; _alg=alg, _maxiters=ode_maxiters)

out = DataFrame(Time_Steps = out_t)
for i in 1:ns
    out[!, Symbol(non_gas_species[i])] = sol_arr[i, :]
end
if use_steady_state && ss_converged
    CSV.write("Profile_ss.csv", out)
else
    CSV.write("Profile_ts.csv", out)
end

try
    profplt = plot(dpi=300, size=(900, 600))
    for i in 1:ns
        sp = String(non_gas_species[i])
        plot!(profplt,
              out.Time_Steps,
              out[!, Symbol(sp)];
              linewidth=3,
              label=sp)
    end
    xlabel!("Time (s)")
    ylabel!("Coverage / Concentration")
    title!("Profiles: adsorbed species")
    plot!(profplt; legend=:outerright)
    if use_steady_state && ss_converged
        png(profplt, "Profiles_ss")
    else
        png(profplt, "Profiles_ts")
    end
catch e
    @warn "Failed to create combined profiles plot: $e"
end

# ==== Step 5b: DRC-ready wrapper ====

function has_arrhenius_columns(df)::Bool
    try
        _get_arr_cols(df; which=:fwd)
        _get_arr_cols(df; which=:rev)
        return true
    catch
        return false
    end
end

function build_k_at_T_arrhenius(df, rx_key, reactions, is_reversible_base, T; Ea_unit::AbstractString="kJ/mol", b::Real=0.0)
    valid_rows = Int[]
    for r in 1:nrow(df)
        s_any = df[r, rx_key]
        s = s_any isa AbstractString ? String(s_any) : string(s_any)
        occursin(ANY_ARROW, s) && push!(valid_rows, r)
    end
    length(valid_rows) == length(reactions) || @warn "Row count for reactions differs; proceeding by order."

    Af_key, Eaf_key = _get_arr_cols(df; which=:fwd)
    Ar_key, Ear_key = _get_arr_cols(df; which=:rev)
    K_key = has_col_case(df, ["K","Keq","K_eq","KEquilibrium","Kconst","Keqconst"]) ?
            find_col_case(df, ["K","Keq","K_eq","KEquilibrium","Kconst","Keqconst"]) : nothing

    Ttyp = typeof(T)
    kf  = Vector{Ttyp}()
    Keq = Vector{Ttyp}()

    for (j, r) in enumerate(valid_rows)
        Af  = _to_float(df[r, Af_key])
        Eaf = _to_float(df[r, Eaf_key])
        Ar  = _to_float(df[r, Ar_key])
        Ear = _to_float(df[r, Ear_key])

        kf_r = arrhenius_k(Af, Eaf, T; Ea_unit=Ea_unit, b=b)
        kr_r = arrhenius_k(Ar, Ear, T; Ea_unit=Ea_unit, b=b)

        push!(kf,  kf_r)
        if K_key !== nothing && is_reversible_base[j]
            # If a thermodynamic K was provided and validated at T0, keep it as the K source.
            push!(Keq, Ttyp(normalize_K(df[r, K_key])))
        else
            push!(Keq, kf_r/kr_r)
        end
    end
    return kf, Keq
end

function log_rate_vec_vs_T(T, analysis_steps, prob, u0, r_input;
                           Ea_unit::AbstractString = get(ENV, "DRC_EA_UNIT", "kJ/mol"),
                           b::Real = 0.0)

    k_fwd_raw_T, K_eq_raw_T = build_k_at_T_arrhenius(df, rx_key, reactions_base, is_reversible_base, T; Ea_unit=Ea_unit, b=b)

    @inbounds for j in eachindex(K_eq_raw_T)
        if !is_reversible_base[j]
            K_eq_raw_T[j] = Inf
        end
    end

    logK_T = map(K -> isfinite(K) ? log(K) : Inf, K_eq_raw_T)
    log_kf_base_T = log.(k_fwd_raw_T)

    p_new = pack_p_mixed(log_kf_base_T, logK_T)
    _prob = remake(prob; p = p_new)

    sol, _ = predict_neuralode(_prob, u0, analysis_steps; _alg=alg, _maxiters=ode_maxiters)
    rate = similar(analysis_steps, eltype(log_kf_base_T))
    @inbounds for i in eachindex(analysis_steps)
        rate_all, _ = compute_reaction_rates(sol[:, i], p_new)
        rate[i] = evaluate_target_reaction(rate_all, r_input)
    end
    return log.(rate)
end

function log_rate_vec_vs_T_noarrhenius(
    T, analysis_steps, prob, u0,
    r_input,
    log_kf_base)

    p_new = pack_p_mixed(log_kf_base, logK_base)
    _prob = remake(prob; p = p_new)

    sol, _ = predict_neuralode(_prob, u0, analysis_steps; _alg=alg, _maxiters=ode_maxiters)
    rate = similar(analysis_steps, eltype(log_kf_base))
    @inbounds for i in eachindex(analysis_steps)
        rate_all, _ = compute_reaction_rates(sol[:, i], p_new)
        rate[i] = evaluate_target_reaction(rate_all, r_input)
    end
    return log.(rate)
end

const _GLOBAL_TARGET = Ref{String}("")

GAS = Set{String}(String(row.Species) for row in eachrow(catalog) if row.Phase == "gas")
logp = Dict(g => log(get(pressure_map, g, 1.0)) for g in GAS)

function pressure_adjustments(reactions, logp, GAS)
    n = length(reactions)
    adj_fwd = zeros(Float64, n)
    adj_rev = zeros(Float64, n)
    @inbounds for j in 1:n
        L, R = reactions[j]
        sF = 0.0
        for (sp, ν) in L
            if sp in GAS; sF += ν * get(logp, sp, 0.0); end
        end
        adj_fwd[j] = sF
        sR = 0.0
        for (sp, ν) in R
            if sp in GAS; sR += ν * get(logp, sp, 0.0); end
        end
        adj_rev[j] = sR
    end
    return adj_fwd, adj_rev
end

function pressure_adjustments_single(reactions, X::AbstractString)
    n = length(reactions)
    aF = zeros(Float64, n)
    aR = zeros(Float64, n)
    @inbounds for j in 1:n
        L, R = reactions[j]
        aF[j] = get(L, String(X), 0.0)
        aR[j] = get(R, String(X), 0.0)
    end
    return aF, aR
end

adj_fwd_base, adj_rev_base = pressure_adjustments(reactions_base, logp, GAS)

logK_base = map(K -> isfinite(K) ? log(K) : Inf, K_eq_raw)

n_rx_total = length(reactions)
@assert length(p) >= n_rx_total
p_tail = p[n_rx_total+1:end]

@inline function pack_p_mixed(log_kf_base::AbstractVector, logK_current::AbstractVector)
    T = eltype(log_kf_base)

    lnkf_eff = log_kf_base .+ T.(adj_fwd_base)
    lnkr_eff = (log_kf_base .- T.(logK_current)) .+ T.(adj_rev_base)

    log_rates = T[]
    @inbounds for j in eachindex(reactions_base)
        push!(log_rates, lnkf_eff[j])
        if is_reversible_base[j]
            push!(log_rates, lnkr_eff[j])
        end
    end

    @assert length(log_rates) == n_rx_total
    return vcat(log_rates, T.(p_tail))
end

# DRC_onesided: independent perturbation on every elementary direction rate constant
# (forward + reverse where reversible), with all other rates fixed.
logk_all_base_fixed = Float64[]
logk_adj_expanded_fixed = Float64[]
active_drc_onesided_idx = Int[]
logk_all_base_raw = Float64[]
logk_adj_expanded_base = Float64[]
drc_onesided_labels = String[]
let expanded_idx = 0
    @inbounds for j in eachindex(reactions_base)
        expanded_idx += 1
        push!(logk_all_base_fixed, log(k_fwd_raw[j]))
        push!(logk_adj_expanded_fixed, adj_fwd_base[j])
        push!(active_drc_onesided_idx, expanded_idx)
        push!(logk_all_base_raw, log(k_fwd_raw[j]))
        push!(logk_adj_expanded_base, adj_fwd_base[j])
        push!(drc_onesided_labels, string(j, "f"))
        if is_reversible_base[j]
            expanded_idx += 1
            logkr = k_rev_raw[j] > 0.0 ? log(k_rev_raw[j]) : -Inf
            push!(logk_all_base_fixed, logkr)
            push!(logk_adj_expanded_fixed, adj_rev_base[j])
            if k_rev_raw[j] > 0.0
                push!(active_drc_onesided_idx, expanded_idx)
                push!(logk_all_base_raw, logkr)
                push!(logk_adj_expanded_base, adj_rev_base[j])
                push!(drc_onesided_labels, string(j, "r"))
            end
        end
    end
end
@assert length(logk_all_base_fixed) == n_rx_total
@assert length(logk_adj_expanded_fixed) == n_rx_total
@assert length(active_drc_onesided_idx) == length(logk_all_base_raw) == length(logk_adj_expanded_base) == length(drc_onesided_labels)

@inline function pack_p_allrates_drc_onesided(logk_all_base::AbstractVector, logk_adj_expanded::AbstractVector)
    T = promote_type(eltype(logk_all_base), eltype(logk_adj_expanded))
    @assert length(logk_all_base) == length(active_drc_onesided_idx)
    @assert length(logk_adj_expanded) == length(active_drc_onesided_idx)
    log_rates_raw = T.(logk_all_base_fixed)
    @inbounds for i in eachindex(active_drc_onesided_idx)
        log_rates_raw[active_drc_onesided_idx[i]] = logk_all_base[i]
    end
    log_rates = log_rates_raw .+ T.(logk_adj_expanded_fixed)
    return vcat(log_rates, T.(p_tail))
end

canonicalize_species(s::AbstractString) = begin
    t = String(s)
    t = replace(t, '\u00D7' => "*")
    t = replace(t, '\u2217' => "*")
    t = replace(t, '\uFE0F' => "")
    t = replace(t, r"\s+" => " ")
    t = strip(t)
    lt = lowercase(t)
    (t == "*" || lt in ("*", "s", "site", "free_site", "freesite")) ? "*" : t
end

"""
    evaluate_target_reaction(rate_all, r_input) -> Real

Compute the net formation rate of `r_input` using the *direction-expanded*
reaction list now stored in `reactions`.

Assumes:
- `reactions` is Vector{Tuple{Dict{String,Float64}, Dict{String,Float64}}} (L,R)
  where reversible steps were expanded into (L,R) and (R,L)
- `rate_all` has length == length(reactions)
"""
function evaluate_target_reaction(rate_all, r_input)
    @assert @isdefined(reactions) "evaluate_target_reaction needs `reactions` in scope"
    n = length(reactions)
    @assert length(rate_all) == n "rate_all must match the expanded reactions length"

    target = canonicalize_species(String(r_input))

    net = zero(eltype(rate_all))
    @inbounds for j in 1:n
        L, R = reactions[j]
        Δν = get(R, target, 0.0) - get(L, target, 0.0)
        if Δν != 0.0
            net += Δν * rate_all[j]
        end
    end
    return net
end

@inline function _pack_p_with_pressures(log_kf_base, adj_fwd, adj_rev, logK, p_tail)
    T = eltype(log_kf_base)
    lnkf_press = log_kf_base .+ adj_fwd
    lnkr_press = (log_kf_base .- logK) .+ adj_rev

    p_new = Vector{T}(undef, n_bidir + length(p_tail))
    @inbounds begin
        p_new[1:n_forward] .= lnkf_press
        p_new[n_forward+1:n_bidir] .= lnkr_press
        p_new[n_bidir+1:end] .= T.(p_tail)
    end
    return p_new
end

valid_species = String.(catalog.Species)

function _pick_target_species(valid::Vector{String})
    target = length(ARGS) >= 2 ? String(ARGS[2]) : get(ENV, "DRC_TARGET", "")

    if isempty(target)
        try
            println("\nEnter target species for analysis (choices: " * join(valid, ", ") * ")")
            print("Target species: "); flush(stdout)
            line = readline(stdin)
            target = strip(line)
            if isempty(target)
                target = (_DEFAULT_TARGET in valid) ? _DEFAULT_TARGET : first(valid)
                println(">> Empty input. Defaulting to $target")
            end
        catch
            target = (_DEFAULT_TARGET in valid) ? _DEFAULT_TARGET : first(valid)
            println(">> Input not available. Defaulting to $target")
        end
    end

    vlower = lowercase.(valid)
    idx = findfirst(==(lowercase(target)), vlower)
    if idx !== nothing
        return valid[idx]
    else
        error("Target species '$target' not found. Options: " * join(valid, ", "))
    end
end

function _read_analysis_grid(default0, default1, defaultn)
    println("\n--- Analysis time grid ---")
    t0   = _read_float("analysis_start", default0)
    t1   = _read_float("analysis_end",   default1)
    npts = _read_int("analysis_points",  defaultn)
    npts < 2 && error("analysis_points must be ≥ 2")
    println("  analysis grid = $(t0) → $(t1)   (n=$(npts))")
    return collect(range(t0, t1; length=npts))
end

function _parse_enthalpy_reference_spec(spec::AbstractString)
    refs = Dict{String,Float64}()
    for raw in split(String(spec), ',')
        item = strip(raw)
        isempty(item) && continue
        parts = occursin("=", item) ? split(item, "="; limit=2) : split(item, ":"; limit=2)
        length(parts) == 2 || error("Malformed enthalpy reference '$item'. Use entries like *=0 or O2:0.")
        refs[strip(parts[1])] = parse(Float64, replace(strip(parts[2]), 'd' => 'e'))
    end
    return refs
end

function read_enthalpy_reference_states(species::Vector{String})
    spec = strip(get(ENV, "DRC_ENTHALPY_REFERENCES", ""))
    if isempty(spec)
        default_spec = "*" in species ? "*=0" : "$(first(species))=0"
        spec = _read_str("Reference enthalpies for state decomposition, e.g. *=0,O2=0,C3H6=0", default_spec)
        isempty(strip(spec)) && (spec = default_spec)
    end
    refs = _parse_enthalpy_reference_spec(spec)
    isempty(refs) && error("At least one enthalpy reference state is required.")
    species_set = Set(species)
    for sp in keys(refs)
        sp in species_set || error("Reference species '$sp' is not in the stable species list: $(join(species, ", "))")
    end
    ref_df = DataFrame(State = collect(keys(refs)), Enthalpy = collect(values(refs)))
    pretty_table(ref_df;
        header = String.(names(ref_df)),
        show_subheader = false,
        tf = tf_unicode_rounded,
        title = "Reference enthalpies for state decomposition",
        alignment = :l,
        crop = :none)
    return refs
end

function side_enthalpy(side::Dict{String,Float64}, H_state::Dict{String,Float64})
    H = 0.0
    for (sp, ν) in side
        H += ν * H_state[sp]
    end
    return H
end

function reconstruct_state_enthalpies(
    reactions_base,
    is_reversible_base,
    stable_species::Vector{String},
    Ea_f::Vector{Float64},
    Ea_r::Vector{Float64},
    reference_H::Dict{String,Float64};
    consistency_tol::Float64 = 1e-3
)
    species_index = Dict(sp => i for (i, sp) in enumerate(stable_species))
    rows = Vector{Vector{Float64}}()
    rhs = Float64[]
    rx_ids = Int[]
    dh_barrier = Float64[]

    for (j, (L, R)) in enumerate(reactions_base)
        is_reversible_base[j] || continue
        row = zeros(Float64, length(stable_species))
        for (sp, ν) in R
            row[species_index[sp]] += ν
        end
        for (sp, ν) in L
            row[species_index[sp]] -= ν
        end
        ΔH = Ea_f[j] - Ea_r[j]
        push!(rows, row)
        push!(rhs, ΔH)
        push!(rx_ids, j)
        push!(dh_barrier, ΔH)
    end

    n_rx_constraints = length(rows)
    for (sp, Href) in reference_H
        row = zeros(Float64, length(stable_species))
        row[species_index[sp]] = 1.0
        push!(rows, row)
        push!(rhs, Href)
    end

    A = reduce(vcat, permutedims.(rows))
    b = rhs
    rnk = rank(A)
    if rnk < length(stable_species)
        @warn "State enthalpy system is underdetermined: rank $(rnk) < $(length(stable_species)). Continuing with least-squares/minimum-norm reconstruction."
    end
    cnd = cond(A)
    if !isfinite(cnd) || cnd > 1e10
        @warn "State enthalpy system is ill-conditioned: cond(A)=$(cnd). Continuing with least-squares reconstruction."
    end

    H_vec = A \ b
    residual = A * H_vec - b
    max_abs_residual = isempty(residual) ? 0.0 : maximum(abs.(residual))
    if max_abs_residual > consistency_tol
        @warn "Thermodynamic inconsistency detected in barrier-derived enthalpies. State-enthalpy decomposition is a directed least-squares reconstruction, not a strictly consistent thermodynamic diagram."
    end

    H_state = Dict(stable_species[i] => H_vec[i] for i in eachindex(stable_species))
    dh_states = Float64[]
    rx_residuals = Float64[]
    for row_i in 1:n_rx_constraints
        j = rx_ids[row_i]
        L, R = reactions_base[j]
        ΔH_states = side_enthalpy(R, H_state) - side_enthalpy(L, H_state)
        push!(dh_states, ΔH_states)
        push!(rx_residuals, ΔH_states - dh_barrier[row_i])
    end
    residual_df = DataFrame(
        Reaction = rx_ids,
        ΔH_from_barriers = dh_barrier,
        ΔH_from_reconstructed_states = dh_states,
        Residual = rx_residuals,
    )

    return H_state, residual_df, max_abs_residual
end

function compute_transition_state_enthalpies(reactions_base, Ea_f, H_state)
    H_TS = zeros(Float64, length(reactions_base))
    for (j, (L, _R)) in enumerate(reactions_base)
        H_TS[j] = side_enthalpy(L, H_state) + Ea_f[j]
    end
    return H_TS
end

function compute_state_drcs(
    drc_onesided,
    drc_onesided_labels,
    reactions_base,
    stable_species::Vector{String},
    analysis_steps
)
    n_t = length(analysis_steps)
    n_rx = length(reactions_base)
    fwd_cols = Dict{Int,Int}()
    rev_cols = Dict{Int,Int}()
    for (col_i, lbl) in enumerate(drc_onesided_labels)
        m = match(r"^(\d+)([fr])$", String(lbl))
        m === nothing && continue
        base_i = parse(Int, m.captures[1])
        if m.captures[2] == "f"
            fwd_cols[base_i] = col_i
        else
            rev_cols[base_i] = col_i
        end
    end

    X_f = zeros(Float64, n_t, n_rx)
    X_r = zeros(Float64, n_t, n_rx)
    for j in 1:n_rx
        haskey(fwd_cols, j) && (X_f[:, j] .= drc_onesided[:, fwd_cols[j]])
        haskey(rev_cols, j) && (X_r[:, j] .= drc_onesided[:, rev_cols[j]])
    end

    X_TS = X_f .+ X_r
    X_stable = zeros(Float64, n_t, length(stable_species))
    for (s_i, sp) in enumerate(stable_species)
        for (j, (L, R)) in enumerate(reactions_base)
            X_stable[:, s_i] .-= get(L, sp, 0.0) .* X_f[:, j]
            X_stable[:, s_i] .-= get(R, sp, 0.0) .* X_r[:, j]
        end
    end

    state_drc_stable = DataFrame(Time = analysis_steps)
    for (s_i, sp) in enumerate(stable_species)
        state_drc_stable[!, Symbol(sp)] = X_stable[:, s_i]
    end

    state_drc_ts = DataFrame(Time = analysis_steps)
    for j in 1:n_rx
        state_drc_ts[!, Symbol("TS$j")] = X_TS[:, j]
    end

    return state_drc_stable, state_drc_ts, X_stable, X_TS
end

while true
    r_input = _pick_target_species(valid_species)
    _GLOBAL_TARGET[] = String(r_input)
    println(">> Target species: ", r_input)

    analysis_steps = _read_analysis_grid(tspan[1], tspan[2], 500)
    println("Running analyses: DRC_Campbell, DRC_onesided (forward/reverse), E_app (if Arrhenius data), and n_X for all gas species.")

    rate_wrapper_drc_campbell_logk = function(log_kf_base)
        p_new = pack_p_mixed(log_kf_base, logK_base)
        _prob = remake(prob; p = p_new)
        times_eval = analysis_steps
        sol, _ = predict_neuralode(_prob, u0, times_eval; _alg=alg, _maxiters=ode_maxiters)
        rate = similar(times_eval, eltype(log_kf_base))
        @inbounds for i in eachindex(times_eval)
            rate_all, _ = compute_reaction_rates(sol[:, i], p_new)
            rate[i] = evaluate_target_reaction(rate_all, r_input)
        end
        return log.(rate)
    end

    rate_wrapper_drc_onesided_logk = function(logk_all_base)
        p_new = pack_p_allrates_drc_onesided(logk_all_base, logk_adj_expanded_base)
        _prob = remake(prob; p = p_new)
        times_eval = analysis_steps
        sol, _ = predict_neuralode(_prob, u0, times_eval; _alg=alg, _maxiters=ode_maxiters)
        rate = similar(times_eval, eltype(logk_all_base))
        @inbounds for i in eachindex(times_eval)
            rate_all, _ = compute_reaction_rates(sol[:, i], p_new)
            rate[i] = evaluate_target_reaction(rate_all, r_input)
        end
        return log.(rate)
    end

    drc_campbell = ForwardDiff.jacobian(rate_wrapper_drc_campbell_logk, log.(k_fwd_raw))
    drc_onesided = ForwardDiff.jacobian(rate_wrapper_drc_onesided_logk, logk_all_base_raw)

    drc_palette = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
    ]
    _base_idx_from_label(lbl::AbstractString) = parse(Int, match(r"^(\d+)", String(lbl)).captures[1])
    _color_for_idx(i::Int) = drc_palette[mod1(i, length(drc_palette))]
    function _subscript_digits(i::Int)
        dmap = Dict(
            '0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
            '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉'
        )
        s = string(i)
        return String([dmap[c] for c in s])
    end
    function _subscript_text(s::AbstractString)
        dmap = Dict(
            '0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
            '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉',
            'x' => 'ₓ', 'X' => 'ₓ'
        )
        out = Char[]
        for c in String(s)
            push!(out, get(dmap, c, c))
        end
        return String(out)
    end

    drc_plot_kwargs = (
        dpi=800,
        size=(1200, 800),
        linewidth=2.2,
        xlabel="Time (s)",
        ylabel="DRC",
        legend=:right,
        tickfontsize=12,
        guidefontsize=14,
        legendfontsize=11,
        framestyle=:box,
        margin=6Plots.mm,
        top_margin=14Plots.mm,
        right_margin=12Plots.mm,
    )

    onesided_fwd_ids = Int[]
    onesided_fwd_cols = Int[]
    onesided_rev_ids = Int[]
    onesided_rev_cols = Int[]
    for (j, lbl) in enumerate(drc_onesided_labels)
        base_i = _base_idx_from_label(lbl)
        if endswith(String(lbl), "f")
            push!(onesided_fwd_ids, base_i)
            push!(onesided_fwd_cols, j)
        elseif endswith(String(lbl), "r")
            push!(onesided_rev_ids, base_i)
            push!(onesided_rev_cols, j)
        end
    end

    plt_drc_campbell = plot(; drc_plot_kwargs...)
    for i in 1:length(reactions_base)
        plot!(plt_drc_campbell, analysis_steps, drc_campbell[:, i];
            color=_color_for_idx(i), label="DRC" * _subscript_digits(i))
    end

    plt_drc_onesided_forward = plot(; drc_plot_kwargs...)
    for (base_i, col_i) in zip(onesided_fwd_ids, onesided_fwd_cols)
        plot!(plt_drc_onesided_forward, analysis_steps, drc_onesided[:, col_i];
            color=_color_for_idx(base_i), label="DRC" * _subscript_digits(base_i))
    end

    plt_drc_onesided_reverse = plot(; drc_plot_kwargs...)
    for (base_i, col_i) in zip(onesided_rev_ids, onesided_rev_cols)
        plot!(plt_drc_onesided_reverse, analysis_steps, drc_onesided[:, col_i];
            color=_color_for_idx(base_i), label="DRC" * _subscript_digits(base_i))
    end

    safe_name = replace(String(r_input), r"[^A-Za-z0-9_]+" => "_")
    png(plt_drc_campbell, "DRC_Campbell_" * safe_name)
    png(plt_drc_onesided_forward, "DRC_onesided_forward_" * safe_name)
    png(plt_drc_onesided_reverse, "DRC_onesided_reverse_" * safe_name)

    drc_campbell_out = DataFrame(Time = analysis_steps)
    drc_onesided_forward_out = DataFrame(Time = analysis_steps)
    drc_onesided_reverse_out = DataFrame(Time = analysis_steps)
    for i in 1:length(reactions_base)
        drc_campbell_out[!, Symbol("DRC_Campbell_$i")] = drc_campbell[:, i]
    end
    for (base_i, col_i) in zip(onesided_fwd_ids, onesided_fwd_cols)
        drc_onesided_forward_out[!, Symbol("DRC_onesided_forward_$(base_i)")] = drc_onesided[:, col_i]
    end
    for (base_i, col_i) in zip(onesided_rev_ids, onesided_rev_cols)
        drc_onesided_reverse_out[!, Symbol("DRC_onesided_reverse_$(base_i)")] = drc_onesided[:, col_i]
    end
    CSV.write("DRC_Campbell_" * safe_name * ".csv", drc_campbell_out)
    CSV.write("DRC_onesided_forward_" * safe_name * ".csv", drc_onesided_forward_out)
    CSV.write("DRC_onesided_reverse_" * safe_name * ".csv", drc_onesided_reverse_out)

    use_arrhenius = has_arrhenius_columns(df)
    if use_arrhenius
        if _LAST_T_ARR[] === nothing
            _LAST_T_ARR[] = _read_float("Base temperature T0 (K)", 500.0)
        end
        if _LAST_EA_UNIT[] === nothing
            _LAST_EA_UNIT[] = _read_choice("Ea unit in CSV (kJ/mol or J/mol)", "kJ/mol"; envkey="DRC_EA_UNIT")
        end
        if _LAST_B_ARR[] === nothing
            _LAST_B_ARR[] = _read_float("Arrhenius exponent b", 0.0; envkey="DRC_B_ARR")
        end

        T0 = _LAST_T_ARR[]
        Ea_unit = _LAST_EA_UNIT[]
        b = _LAST_B_ARR[]

        invT0 = 1.0 / T0

        # Direct variable transform: x = 1/T.
        # E_app = -R * d(ln r)/d(1/T)
        f_invT = function (xvec)
            x = xvec[1]
            T = 1.0 / x
            return log_rate_vec_vs_T(
                T, analysis_steps, prob, u0, r_input;
                Ea_unit = Ea_unit,
                b = b)
        end

        J_invT = ForwardDiff.jacobian(f_invT, [invT0])
        dlnr_dinvT = vec(J_invT)
        
        Eapp = (-R_gas) .* dlnr_dinvT

        out_unit = lowercase(String(Ea_unit)) in ("kj/mol","kjmol","kj") ? "kJ/mol" : "J/mol"
        scale = (lowercase(out_unit) in ("kj/mol","kj")) ? 1e-3 : 1.0
        Eapp_disp = scale .* Eapp

        safe_target = replace(String(r_input), r"[^A-Za-z0-9_]+" => "_")

        _, Eaf_key = _get_arr_cols(df; which=:fwd)
        _, Ear_key = _get_arr_cols(df; which=:rev)
        Ea_f_raw = [_to_float(df[r, Eaf_key]) for r in 1:length(reactions_base)]
        Ea_r_raw = [_to_float(df[r, Ear_key]) for r in 1:length(reactions_base)]
        # Ea values loaded from CSV are already in Ea_unit; use them directly for weighting.
        # (E_app itself is converted to out_unit via `scale` above.)
        Ea_f_disp = Ea_f_raw
        Ea_r_disp = Ea_r_raw

        Eapp_onesided_forward_weighted = zeros(length(analysis_steps))
        for (base_i, col_i) in zip(onesided_fwd_ids, onesided_fwd_cols)
            Eapp_onesided_forward_weighted .+= drc_onesided[:, col_i] .* Ea_f_disp[base_i]
        end

        Eapp_onesided_reverse_weighted = zeros(length(analysis_steps))
        for (base_i, col_i) in zip(onesided_rev_ids, onesided_rev_cols)
            Eapp_onesided_reverse_weighted .+= drc_onesided[:, col_i] .* Ea_r_disp[base_i]
        end
        Eapp_onesided_weighted = Eapp_onesided_forward_weighted .+ Eapp_onesided_reverse_weighted

        state_decomp_ok = false
        Eapp_state_weighted = fill(NaN, length(analysis_steps))
        max_state_enthalpy_residual = NaN

        # Correct Campbell-style E_app decomposition: stable-state enthalpies plus
        # transition-state enthalpies weighted by state DRCs. The old direct
        # Campbell DRC × forward barrier curve is intentionally not used here.
        if !any(is_reversible_base)
            println(">> Skipping state-enthalpy decomposition: no reversible reactions available for enthalpy reconstruction.")
        else
            stable_species = String.(display_species)
            reference_H = read_enthalpy_reference_states(stable_species)
            try
                residual_tol = lowercase(out_unit) == "kj/mol" ? 1e-3 : 1.0
                H_state, enthalpy_residual_df, max_state_enthalpy_residual = reconstruct_state_enthalpies(
                    reactions_base,
                    is_reversible_base,
                    stable_species,
                    Float64.(Ea_f_disp),
                    Float64.(Ea_r_disp),
                    reference_H;
                    consistency_tol = residual_tol)
                H_TS = compute_transition_state_enthalpies(reactions_base, Ea_f_disp, H_state)
                state_drc_stable, state_drc_ts, X_stable, X_TS = compute_state_drcs(
                    drc_onesided,
                    drc_onesided_labels,
                    reactions_base,
                    stable_species,
                    analysis_steps)

                Eapp_state_weighted .= 0.0
                for (s_i, sp) in enumerate(stable_species)
                    Eapp_state_weighted .+= X_stable[:, s_i] .* H_state[sp]
                end
                for j in eachindex(reactions_base)
                    Eapp_state_weighted .+= X_TS[:, j] .* H_TS[j]
                end

                state_enthalpy_out = DataFrame(State = String[], Type = String[], Enthalpy = Float64[])
                for sp in stable_species
                    push!(state_enthalpy_out, (sp, "stable", H_state[sp]))
                end
                for j in eachindex(H_TS)
                    push!(state_enthalpy_out, ("TS$j", "transition_state", H_TS[j]))
                end
                CSV.write("State_enthalpies_" * safe_target * ".csv", state_enthalpy_out)
                CSV.write("State_enthalpy_consistency_" * safe_target * ".csv", enthalpy_residual_df)
                CSV.write("State_DRC_stable_" * safe_target * ".csv", state_drc_stable)
                CSV.write("State_DRC_TS_" * safe_target * ".csv", state_drc_ts)

                summary_df = DataFrame(
                    Quantity = [
                        "Workflow E_app at final analysis time",
                        "One-sided barrier weighted E_app at final analysis time",
                        "State-enthalpy weighted E_app at final analysis time",
                        "State - Workflow difference",
                        "One-sided - Workflow difference",
                        "Max enthalpy consistency residual",
                    ],
                    Value = [
                        Eapp_disp[end],
                        Eapp_onesided_weighted[end],
                        Eapp_state_weighted[end],
                        Eapp_state_weighted[end] - Eapp_disp[end],
                        Eapp_onesided_weighted[end] - Eapp_disp[end],
                        max_state_enthalpy_residual,
                    ],
                )
                pretty_table(summary_df;
                    header = String.(names(summary_df)),
                    show_subheader = false,
                    tf = tf_unicode_rounded,
                    title = "State enthalpy decomposition summary",
                    alignment = :l,
                    crop = :none)
                state_decomp_ok = true
            catch e
                @warn "Skipping state-enthalpy decomposition after reconstruction failure: $e"
            end
        end

        eapp_time_span = maximum(analysis_steps) - minimum(analysis_steps)
        eapp_dx = eapp_time_span > 0 ? 0.006 * eapp_time_span : 0.0
        eapp_onesided_jitter = [0.75 + 0.25 * sin(2 * pi * (i - 1) / 5) for i in eachindex(analysis_steps)]
        eapp_state_jitter = [0.75 + 0.25 * cos(2 * pi * (i - 1) / 7) for i in eachindex(analysis_steps)]
        eapp_time_onesided = collect(analysis_steps) .- eapp_dx .* eapp_onesided_jitter
        eapp_time_state = collect(analysis_steps) .+ eapp_dx .* eapp_state_jitter

        plt_eapp = plot(; dpi=800, size=(1200, 800),
            xlabel="Time (s)",
            ylabel="E_app ($out_unit)",
            legend=:outerright,
            tickfontsize=12,
            guidefontsize=14,
            legendfontsize=11,
            framestyle=:box,
            title="Apparent activation energy of $(String(r_input))")

        plot!(plt_eapp, analysis_steps, Eapp_disp;
            linewidth=2.8,
            color=:black,
            label="Workflow E_app")
        scatter!(plt_eapp, eapp_time_onesided, Eapp_onesided_weighted;
            label="One-sided barrier weighted",
            color="#2ca02c",
            marker=:circle,
            markersize=4,
            markerstrokewidth=0.8,
            alpha=0.9)
        if state_decomp_ok
            scatter!(plt_eapp, eapp_time_state, Eapp_state_weighted;
                label="State-enthalpy weighted",
                color="#1f77b4",
                marker=:diamond,
                markersize=4,
                markerstrokewidth=0.8,
                alpha=0.9)
        end

        png(plt_eapp, "Eapp_" * safe_target)

        Eapp_out = DataFrame(
            Time = analysis_steps,
            E_app_workflow = Eapp_disp,
            E_app_onesided_barrier_weighted = Eapp_onesided_weighted,
            E_app_state_enthalpy_weighted = Eapp_state_weighted,
            E_app_state_minus_workflow = Eapp_state_weighted .- Eapp_disp,
            E_app_onesided_minus_workflow = Eapp_onesided_weighted .- Eapp_disp,
            max_state_enthalpy_residual = fill(max_state_enthalpy_residual, length(analysis_steps)),
        )
        CSV.write("Eapp_" * safe_target * ".csv", Eapp_out)
    else
        println(">> Skipping E_app: Arrhenius columns not detected in input CSV.")
    end

    gas_species_all = collect(GAS)
    if isempty(gas_species_all)
        println(">> No gas-phase species found; skipping apparent reaction orders.")
    else
        log_kf_base = log.(k_fwd_raw)

        function pack_p_mixed_custom(log_kf_base::AbstractVector, logK_current::AbstractVector,
                                     adjF::AbstractVector, adjR::AbstractVector, p_tail)

            T = promote_type(eltype(log_kf_base), eltype(logK_current), eltype(adjF), eltype(adjR), eltype(p_tail))

            log_kfT = T.(log_kf_base)
            logKT   = T.(logK_current)
            adjFT   = T.(adjF)
            adjRT   = T.(adjR)
            tailT   = T.(p_tail)

            lnkf_eff = log_kfT .+ adjFT
            lnkr_eff = (log_kfT .- logKT) .+ adjRT

            log_rates = T[]
            @inbounds for j in eachindex(reactions_base)
                push!(log_rates, lnkf_eff[j])
                if is_reversible_base[j]
                    push!(log_rates, lnkr_eff[j])
                end
            end

            @assert length(log_rates) == length(reactions)
            return vcat(log_rates, tailT)
        end

        for X in gas_species_all
            try
                aF_X, aR_X = pressure_adjustments_single(reactions_base, X)

                f_logpX = function(svec)
                    s = svec[1]

                    adjF = adj_fwd_base .+ s .* aF_X
                    adjR = adj_rev_base .+ s .* aR_X

                    p_new = pack_p_mixed_custom(log_kf_base, logK_base, adjF, adjR, p_tail)
                    _prob = remake(prob; p = p_new)

                    sol, _ = predict_neuralode(_prob, u0, analysis_steps; _alg=alg, _maxiters=ode_maxiters)
                    rate = Vector{typeof(s)}(undef, length(analysis_steps))
                    @inbounds for i in eachindex(analysis_steps)
                        rate_all, _ = compute_reaction_rates(sol[:, i], p_new)
                        rate[i] = evaluate_target_reaction(rate_all, r_input)
                    end

                    return log.(rate)
                end

                J  = ForwardDiff.jacobian(f_logpX, [0.0])
                nX = vec(J)

                safe_target = replace(String(r_input), r"[^A-Za-z0-9_]+" => "_")
                safe_X      = replace(String(X),       r"[^A-Za-z0-9_]+" => "_")

                nx_label = "n₍" * _subscript_text(String(X)) * "₎"
                plt_nx = plot(analysis_steps, nX;
                        dpi=800, size=(1200,800),
                        linewidth=2.2, color="#1f77b4", label=nx_label,
                        xlabel="Time (s)", ylabel="Apparent order nₓ",
                        legend=:outerright,
                        tickfontsize=12,
                        guidefontsize=14,
                        legendfontsize=11,
                        framestyle=:box,
                        title = "Order of $(X) on r($(safe_target))")
                png(plt_nx, "nX_$(safe_X)_on_$(safe_target)")

                nX_out = DataFrame(Time = analysis_steps, n_X = nX)
                CSV.write("nX_$(safe_X)_on_$(safe_target).csv", nX_out)
            catch e
                println(">> Skipping n_X for gas species $(X): $(typeof(e))")
            end
        end
    end

    resp = _read_str("Run another target? Press ESC then Enter to exit; Enter to continue", "")
    if !isempty(resp)
        any(c -> c == Char(27), resp) && break
        (lowercase(strip(resp)) in ("n","no","q","quit","exit")) && break
    end
end

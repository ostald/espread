
function save_commit_hash(res_dir)
    # saves current commit hash for reproducability
    current_commit = read(`git rev-parse HEAD`, String)
    open(joinpath(res_dir, "git_commit.txt"),"a") do io
        println(io,"current_commit=",current_commit)
    end
end

function auto_commit(res_dir, commit_message)
    run(`git add -u`)
    commit_string = "auto_commit: " * commit_message
    current_commit = read(`git commit -m $commit_string`, String)
    #return current_commit
    print(current_commit)
end

using Serialization

"""
    load_ee_comparison(res_file)

Load a result file written by collison_free_electron.jl (with vs. without e-e energy loss).
Returns a NamedTuple with the header (E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals),
the traces (st_on, r_on, v_on, t_on, st_off, r_off, v_off, t_off; r and v are 3 x n_steps),
and per-step altitude [m] and energy [eV] (alt_on, alt_off, E_on, E_off).
Needs constants.jl (altitude, E_ev) to be included.
"""
function load_ee_comparison(res_file)
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, tr = open(res_file) do io
        (deserialize(io)..., deserialize(io))
    end
    alt_on = altitude.(collect.(eachcol(tr.r_on)))
    alt_off = altitude.(collect.(eachcol(tr.r_off)))
    E_on = E_ev.(norm.(eachcol(tr.v_on)))
    E_off = E_ev.(norm.(eachcol(tr.v_off)))
    return (; E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, tr..., alt_on, alt_off, E_on, E_off)
end
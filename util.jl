
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
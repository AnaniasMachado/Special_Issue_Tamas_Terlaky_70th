
function rename_files(folder::String, old_substring::String, new_substring::String)
    files = readdir(folder)

    for file in files
        old_path = joinpath(folder, file)
        
        if isfile(old_path)
            new_name = replace(file, old_substring => new_substring)
            new_path = joinpath(folder, new_name)

            if new_name != file
                println("Renaming '$file' to '$new_name'")
                mv(old_path, new_path)
            end
        end
    end
end

folder = "./instances/square_dense_test"
old_substring = "experiment_5_matrix"
new_substring = "A"

rename_files(folder, old_substring, new_substring)
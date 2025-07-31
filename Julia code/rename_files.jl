
function rename_files(folder::String, old_substring::String, new_substring::String)
    files = readdir(folder)

    for file in files
        old_path = joinpath(folder, file)
        
        if isfile(old_path)
            new_name = replace(file, old_substring => new_substring)
            new_path = joinpath(folder, new_name)

            if new_name != file
                println("Renomeando '$file' para '$new_name'")
                mv(old_path, new_path)
            end
        end
    end
end

folder = "./Solutions/DRS_Experiment_2/DRS_FP"
old_substring = "PLS"
new_substring = "P13"

rename_files(folder, old_substring, new_substring)
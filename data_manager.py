import os
import shutil


def delete_files_in_subfolders(directory):
    # Walk through the directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            os.remove(file_path)  # Delete the file


def copy_all_to_new_folder(source_directory, target_base_directory, new_folder_name):
    # Construct the full target directory path
    target_directory = os.path.join(target_base_directory, new_folder_name)

    # Ensure the target directory does not already exist to avoid errors
    if os.path.exists(target_directory):
        raise Exception(
            f"The folder '{target_directory}' already exists. Please use a different name or delete the existing folder."
        )
    # Copy the entire directory tree from the source to the target location
    shutil.copytree(source_directory, target_directory)

    return target_directory


def manage_files(
    source_directory,
    target_base_directory,
    new_folder_name,
    file_to_copy="ziji_param.json",
):
    # Copy all files and directories from the source to the new folder
    new_folder_path = copy_all_to_new_folder(
        source_directory, target_base_directory, new_folder_name
    )

    # Copy the specific file to the new folder
    extra_file_path = file_to_copy
    if os.path.exists(extra_file_path):
        shutil.copy2(extra_file_path, new_folder_path)
    else:
        print(f"Warning: The file '{file_to_copy}' was not found and thus not copied.")

    # Delete all files in the target base directory (but not the directories themselves)
    delete_files_in_subfolders(source_directory)


if __name__ == "__main__":

    source_directory_path = "data/data_out/"
    target_base_directory = "ziji_disk"
    new_folder_name = "NewFolder"

    manage_files(source_directory_path, target_base_directory, new_folder_name)

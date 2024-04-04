# Read the file containing file names to keep
$keepFileNames = Get-Content -Path "instance_names.txt"

# Specify the directory to clean up
$directoryPath = ".\instances\"

# Get all files in the directory
$filesInDirectory = Get-ChildItem -Path $directoryPath

# Iterate through each file in the directory
foreach($file in $filesInDirectory){
    # If the file's name is not in the list of names to keep, delete the file
    if($keepFileNames -notcontains ("dantzig42" -or
                                    "swiss42" -or 
                                    "att48" -or
                                    "gr48" -or
                                    "hk48" -or
                                    "eil51 " -or
                                    "berlin52" -or
                                    "brazil58" -or
                                    "st70 " -or
                                    "eil76 " -or
                                    "pr76" -or
                                    "pr76r" -or
                                    "gr96" -or
                                    "rat99" -or
                                    "kroA100" -or
                                    "kroB100" -or
                                    "kroC100" -or
                                    "kroD100" -or
                                    "kroE100" -or
                                    "rd100" -or
                                    "eil101" -or
                                    "lin105x" -or
                                    "pr107")){
        Remove-Item -Path $file.FullName
    }
}
# Data Generation Mat Files Procedure

1. Convert ROS data files (.csv) to a session text file using `csv_to_txt.m`. The session text file will be saved in `/RawData_T1_400`

Ensure the triggers matches the one in ROS triggers file and that the number of triggers is in multiple of 3

2. Generate eyelink and unity mat files to be used for data generation using `unityfile.m` and `eyelink.m`

## UnityFile.m

Run this command to generate `unityfile.mat`
```
unityfile ('auto','save','redo)
```

## Eyelink.m

Run this command to generate `eyelink.mat`
```
eyelink ('auto','save','redo')
```

## Pointers

### Unityfile

- Adjust FileLineOffset in `Line 15` to match the number of lines to offset to the first data point, usually denoted by `Trigger Version 84`

### Eyelink 

- If `Trigger Version 84` is not seen in the eyelink data file (.edf), the variable `s` in `line 197` can be modified to match the start of the first trial. This can be implemented dynamically in the future.
- If the triggers are differently named, modify `completeData.m` in `line 206 - 229` to match the triggers

For both matlab files, do run the above scripts in the `/SessionXX` folder

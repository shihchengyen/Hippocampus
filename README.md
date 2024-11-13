# Hippocampus

## Usage

### Place cell
The `vmpc` object contains information about a cell's encoding of space. To create it, navigate to the cell directory and instantiate the object

```matlab
cd <cell directory>
vmp = vmpc('auto');
```

where `<cell directory>` should be replaced by an actual directory, e.g. "session01/array01/channel001/cell01". To plot a smoothed version of the cell's place field, run

```matlab
plot(vmp, 1, 'MapOnly')
```

This should pop-up a figure contained the smoothed place field of this cell, i.e. a map of the locations in the maze where this cell's spikes predominantly occured.
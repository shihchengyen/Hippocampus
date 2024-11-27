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

### View cell

A view cell encodes information about what the subject is looking at. We compute the view represention of a cell by creating a `vmsv` object.

```matlab
cd <cell directory>
vmv = vmsv('auto','Prefix','1');
```

We can get a 3D representation of the view information contained in this cell by doing

```matlab
plotmap(vmv.data.maps_adsm, 'view')
```
This should pop-up a so-called exploded view, showing the density of spikes elicited when the subject looked at the walls, the pillars, the floor and the ceiling.
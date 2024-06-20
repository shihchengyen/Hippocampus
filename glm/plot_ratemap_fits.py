import os
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Suppress all warnings
import warnings
warnings.filterwarnings('ignore')


# Define key constants here
num_place, num_hd, num_view = 1600, 60, 5122


def cosine_similarity(arr1, arr2):
    return np.nansum(np.multiply(arr1, arr2)) / (np.sqrt(np.nansum(np.multiply(arr1, arr1))) * np.sqrt(np.nansum(np.multiply(arr2, arr2))))   


def model_ratemaps(hc_file, var_type):
    # Load in glm_hardcastle_results.mat file
    hc_results = h5py.File(hc_file).get('hc_results')
    tbin_size = hc_results.get('tbin_size')[0,0]
    num_folds = int(hc_results.get('num_folds')[0,0])
    model_class = int(np.nan_to_num(hc_results.get('classification')[0,0]))
    params_consol = hc_results.get('params_consol')
    similarity_scores = hc_results.get('similarity_scores')[0]
    model_params = list()
    for i in range(params_consol.shape[0]):
        row = list()
        for j in range(params_consol.shape[1]):
            row.append(np.array(hc_results[params_consol[i,j]]).flatten())
        model_params.append(row)
    model_scores = list()
    for i in range(similarity_scores.shape[0]):
        model_scores.append(np.array(hc_results[similarity_scores[i]]))

    # Extract model-fitted parameters based on classification
    if model_class == 1:  # phv
        place_params = np.array(model_params[0])[:,:num_place].T
        hd_params = np.array(model_params[0])[:,num_place:num_place+num_hd].T
        view_params = np.array(model_params[0])[:,num_place+num_hd:].T

        place_scores = model_scores[0][0,:]
        hd_scores = model_scores[0][1,:]
        view_scores = model_scores[0][2,:]

    elif model_class == 2:  # ph
        place_params = np.array(model_params[1])[:,:num_place].T
        hd_params = np.array(model_params[1])[:,num_place:].T
        view_params = np.array(model_params[6]).T

        place_scores = model_scores[1][0,:]
        hd_scores = model_scores[1][1,:]
        view_scores = model_scores[6].flatten()

    elif model_class == 3:  # pv
        place_params = np.array(model_params[2])[:,:num_place].T
        hd_params = np.array(model_params[5]).T
        view_params = np.array(model_params[2])[:,num_place:].T

        place_scores = model_scores[2][0,:]
        hd_scores = model_scores[5].flatten()
        view_scores = model_scores[2][1,:]

    elif model_class == 4:  # hv
        place_params = np.array(model_params[4]).T
        hd_params = np.array(model_params[3])[:,:num_hd].T
        view_params = np.array(model_params[3])[:,num_hd:].T

        place_scores = model_scores[4].flatten()
        hd_scores = model_scores[3][0,:]
        view_scores = model_scores[3][1,:]

    else:  # single-variable models or unclassified
        place_params = np.array(model_params[4]).T
        hd_params = np.array(model_params[5]).T
        view_params = np.array(model_params[6]).T

        place_scores = model_scores[4].flatten()
        hd_scores = model_scores[5].flatten()
        view_scores = model_scores[6].flatten()

    # Get firing rate map from model-fitted parameters
    if var_type == 'place':
        model_params = place_params
        model_scores = place_scores
    elif var_type == 'headdirection':
        model_params = hd_params
        model_scores = hd_scores
    elif var_type == 'spatialview':
        model_params = view_params
        model_scores = view_scores

    fr_model = np.exp(model_params)/tbin_size
    return fr_model, model_scores


def session_raw_ratemap(var_type, num_folds):
    vmobj_filenames = {'place': 'vmpc', 'headdirection': 'vmhd', 'spatialview': 'vmsv'}
    vmobj_varnames = {'place': 'vmp', 'headdirection': 'vmd', 'spatialview': 'vms'}

    # Load in vmobj.mat file
    vmobj_data = h5py.File(f'../{vmobj_filenames[var_type]}.mat').get(vmobj_varnames[var_type]).get('data')
    fr_data = np.array(vmobj_data.get('maps_raw')).flatten()
    bin_filt = ~np.isnan(fr_data)
    fr_data = [fr_data for _ in range(num_folds)]
    return fr_data, bin_filt


def session_smoothed_ratemap(var_type, num_folds):
    vmobj_filenames = {'place': 'vmpc', 'headdirection': 'vmhd', 'spatialview': 'vmsv'}
    vmobj_varnames = {'place': 'vmp', 'headdirection': 'vmd', 'spatialview': 'vms'}

    # Load in vmobj.mat file
    vmobj_data = h5py.File(f'../{vmobj_filenames[var_type]}.mat').get(vmobj_varnames[var_type]).get('data')
    fr_data = np.array(vmobj_data.get('maps_adsm')).flatten()
    bin_filt = ~np.isnan(fr_data)
    fr_data = [fr_data for _ in range(num_folds)]
    return fr_data, bin_filt


def training_raw_ratemap(var_type, num_folds):
    # Load in vmpvData.mat file
    vmpv_data = h5py.File('vmpvData.mat').get('vmpvData')
    bin_stc = np.array(vmpv_data.get('bin_stc')).T
    tbin_size = vmpv_data.get('tbin_size')[0,0]
    place_good_bins = np.array(vmpv_data.get('place_good_bins'), dtype=int).flatten()-1
    view_good_bins = np.array(vmpv_data.get('view_good_bins'), dtype=int).flatten()-1
    
    # Get params based on var_type
    if var_type == 'place':
        num_params = num_place
        pcol = 1
        good_bins = place_good_bins
    elif var_type == 'headdirection':
        num_params = num_hd
        pcol = 2
        good_bins = np.arange(num_params)
    elif var_type == 'spatialview':
        num_params = num_view
        pcol = 3
        good_bins = view_good_bins

    # Generate firing rate maps per bin, using training data from each fold
    fold_edges = np.round(np.linspace(0, bin_stc.shape[0], (5*num_folds)+1))
    fr_data = list()
    for k in range(num_folds):
        # Generate training datapoints used for each fold
        test_ind  = np.hstack([np.arange(fold_edges[k], fold_edges[k+1]), np.arange(fold_edges[k+num_folds], fold_edges[k+num_folds+1]), np.arange(fold_edges[k+2*num_folds], fold_edges[k+2*num_folds+1]),\
                    np.arange(fold_edges[k+3*num_folds], fold_edges[k+3*num_folds+1]), np.arange(fold_edges[k+4*num_folds], fold_edges[k+4*num_folds+1])])
        train_ind = np.setdiff1d(np.arange(bin_stc.shape[0]), test_ind)
        bin_stc_fold = bin_stc[train_ind,:]
        # Generate distribution of spike counts per observation for each bin
        fr_fold = [list() for _ in range(num_params)]
        for row in range(bin_stc_fold.shape[0]):
            bin_num = int(bin_stc_fold[row, pcol])-1
            fr_fold[bin_num].append(bin_stc_fold[row, 4])
        fr_fold = list(map(np.array, fr_fold))
        fr_fold = np.array(list(map(lambda arr: np.sum(arr)/(tbin_size*arr.shape[0]), fr_fold)))
        fr_data.append(fr_fold)
    bin_filt = good_bins

    return fr_data, bin_filt


def plot_ratemaps(fr_model, fr_data, model_scores, var_type, ax_range=[]):
    num_folds = fr_model.shape[1]
    var_names = {'place': 'place', 'headdirection': 'hd', 'spatialview': 'view'}

    # Plot vmobj adsmoothed ratemap and model-fitted ratemap together
    fig, axes = plt.subplots(num_folds, 1, figsize=(25, 10), sharey=True)
    if var_type == 'place':
        line_width = 1.2
    elif var_type == 'spatialview':
        line_width = 0.8
    else:
        line_width = 1
    for i in range(num_folds):
        ax = axes[i]
        ax.plot(np.arange(fr_model.shape[0]), fr_model[:,i], label='glm fit', color='C1', linewidth=line_width)
        ax.bar(np.arange(fr_data[i].shape[0]), fr_data[i], label='fold data', color='C0')
        ax.text(0.995, 0.9, f'Ratemap fit: {model_scores[i]:.6f}', transform=ax.transAxes, ha='right', va='top')
        ax.set_xticks([])
        if len(ax_range) == 2:
            ax.set_ylim(ax_range[0], ax_range[1])
        else:
            ax.set_ylim(0, 1.5*np.max([np.nanpercentile(fr_model, 99), np.nanpercentile(fr_data, 99)]))
        
    axes[0].legend(loc='upper left')
    fig.text(0.5, 0.09, 'Bin #', ha='center')
    fig.text(0.105, 0.5, 'Firing rate (Hz)', va='center', rotation='vertical')
    fig.suptitle(f'{var_names[var_type].capitalize()} firing rate maps per fold', y=0.905)

    plt.savefig(f'{var_names[var_type]}_ratemap_fits.png', bbox_inches='tight')
    #plt.show(block=False)
    return


def main():
    # Get command line args
    if len(sys.argv) >= 3:
        var_type = sys.argv[1]
        reference = sys.argv[2]
        hc_file = sys.argv[3]
    elif len(sys.argv) == 2:
        var_type = sys.argv[1]
        reference = sys.argv[2]
        hc_file = 'glm_hardcastle_results.mat'
    elif len(sys.argv) == 1: 
        var_type = sys.argv[1]
        reference = 'training'
        hc_file = 'glm_hardcastle_results.mat'
    else:
        # Need at least 1 cli argument: var_type
        return
    
    # Other parameters manually specified here:
    recompute = True
    scale = True
    filt = False

    # Load in glm_hardcastle_results.mat file and get firing ratemap from model-fitted parameters
    print(f'Loading {hc_file} and generating model firing rate map for {var_type}...')
    fr_model, model_scores = model_ratemaps(hc_file, var_type)
    num_folds = fr_model.shape[1]

    # Load in reference firing ratemaps (vmobj raw ratemap/vmobj smoothed ratemap/vmpvData training folds)
    print(f'Loading reference firing rate map ({reference})...')
    if reference == 'raw':
        fr_data, bin_filt = session_raw_ratemap(var_type, num_folds)
    elif reference == 'smoothed':
        fr_data, bin_filt = session_smoothed_ratemap(var_type, num_folds)
    elif reference == 'training':
        fr_data, bin_filt = training_raw_ratemap(var_type, num_folds)

    # Whether to use ratemap fit scores from glm_hardcastle_results.mat or recompute with current data
    if recompute:
        # Currently only have cosine similarity implemented
        model_scores = [cosine_similarity(fr_model[:,i], fr_data[i]) for i in range(num_folds)]

    # Whether or not to scale firing rates by the norm of the firing rate vector
    if scale:
        fr_model = fr_model/np.nansum(fr_model)
        fr_data = fr_data/np.nansum(fr_data)

    # Whether or not to filter out unoccupied bins when plotting
    if filt:
        fr_model = fr_model[bin_filt,:]
        fr_data = [fr_fold[bin_filt] for fr_fold in fr_data]

    # Plot vmobj adsmoothed ratemap and model-fitted ratemap together
    print(f'Generating and saving plot...')
    plot_ratemaps(fr_model, fr_data, model_scores, var_type)
    return


if __name__ == '__main__':
    main()

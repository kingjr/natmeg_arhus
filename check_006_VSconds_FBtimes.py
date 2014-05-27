
# coding: utf-8

#### import settings


import numpy as np
import mne

import os
# directories ISIS
data_path = "/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/" + \
                "scratch/tsss_initial/006_HEN/"
eve_path = "/projects/MINDLAB2013_01-MEG-AttentionEmotionVisualTracking/" + \
                "scratch/events.fif/006_HEN/raw/"
                
data_path = "/media/VSCOND_MEG/Arhus_data/tSSS/006_HEN/"
eve_path = "/media/VSCOND_MEG/Arhus_data/events.fif/006_HEN/raw/"
            
os.chdir(data_path)

n_jobs = 4 # number of processors to use

# epoch variables
tmin, tmax = -0.2, 0.350
baseline = (-0.2, -0.050)  # baseline time
#reject = dict(mag=4e-12, grad=4000e-13)

event_ids = {"FB_neu3": 10, "FB_ang3": 11, "FB_neu10": 20, "FB_ang10": 21}  

#### Processing raw files


for cond in ["VS_1b_1","VS_1b_2"]:
    #raw = mne.fiff.Raw(data_path +
    #                          "VS_1b_1_tsss_mc.fif", preload=False)
    raw = mne.fiff.Raw(data_path + cond + '_tsss_mc.fif', preload=False) 
    #                          "VS_1b_2_tsss_mc.fif", preload=False)
    
    picks_ana = mne.fiff.pick_channels(raw.info['ch_names'], include='MISC001')
    picks = mne.fiff.pick_types(raw.info, meg='grad', eeg=False, eog=True,
                                       emg=True, ecg=True, misc=True, exclude='bads')
    
    
    events_uc = mne.find_events(raw, stim_channel='STI101', min_duration=0.002)
    epochs_uc = mne.Epochs(raw, events_uc, event_ids, tmin, tmax,
                               picks=picks, baseline=baseline,
                               preload=True, reject=None)
    epochs_uc_ana = mne.Epochs(raw, events_uc, event_ids, tmin, tmax,
                               picks=picks_ana, baseline=baseline,
                               preload=True, reject=None)
    
    evoked_uc_FB = epochs_uc[['FB_neu3','FB_ang3','FB_neu10','FB_ang10']].average()
    evoked_uc_FB_ana = epochs_uc_ana[['FB_neu3','FB_ang3','FB_neu10','FB_ang10']].average(picks=[0])
    
    events_c = mne.read_events(eve_path + cond + '-eve.fif') #'"VS_1b_2-eve.fif")
    epochs_c = mne.Epochs(raw, events_c, event_ids, tmin, tmax,
                               picks=picks, baseline=baseline,
                               preload=True, reject=None)
    epochs_c_ana = mne.Epochs(raw, events_c, event_ids, tmin, tmax,
                               picks=picks_ana, baseline=baseline,
                               preload=True, reject=None)
    
    evoked_c_FB = epochs_c[['FB_neu3','FB_ang3','FB_neu10','FB_ang10']].average()
    evoked_c_FB_ana = epochs_c_ana[['FB_neu3','FB_ang3','FB_neu10','FB_ang10']].average(picks=[0])

    mne.viz.plot_image_epochs(epochs_uc_ana, picks=[0])
    mne.viz.plot_image_epochs(epochs_c_ana, picks=[0])
    
    mne.viz.plot_topo([evoked_uc_FB, evoked_c_FB],
                      color=["green", "yellow"], title=cond + ":UnCorr vs. Corr")

    
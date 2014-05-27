import mne
import scipy.io as sio
import sys as sys
file_events = str(sys.argv[2])
file_data = str(sys.argv[1])
events = mne.read_events(file_events)
raw = mne.fiff.Raw(file_data)
start = raw.first_samp
print "start:"+str(start)
print "events[0]:"+str(events[0,0])+"->"+str(events[0,0]-start)

events[:,0] = events[:,0]-start
sio.savemat(file_events[0:-3] + "mat", {'events':events})
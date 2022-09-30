# Module for importing Intan data saved in "One File Per Signal Type" format
import spikeinterface.extractors as se
import spikeinterface.core as sc

def read_intan_dat(data_path):

    info_file = data_path/'info.rhd'
    rec_info = se.read_intan(info_file, stream_id='0')

    sampling_frequency = rec_info.get_sampling_frequency()
    num_channels = rec_info.get_num_channels()
    # num_segments = rec_info.get_num_segments() # Should be 1 - could add test?

    # Loading the .dat as a generic binary
    recording_file  = data_path/'amplifier.dat'
    info = rec_info.to_dict()['properties']
    gain_to_uV = info['gain_to_uV'][0]  # Should be 0.195
    offset_to_uV = 0                    # Don't use info['offset_to_uV']! This should be 0
    dtype = "int16"
    time_axis = 0
    recording = sc.read_binary(recording_file, num_chan=num_channels, sampling_frequency=sampling_frequency,
                            dtype=dtype, gain_to_uV=gain_to_uV, offset_to_uV=offset_to_uV, 
                            time_axis=time_axis)
    # recording.annotate(is_filtered=False) # This is the default
    return recording
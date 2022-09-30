import numpy as np
import probeinterface as pi

# Cambridge neurotech probe 'ASSY-236-H6'
# If given info_rhd_file, will automatically remove channels that weren't recorded!
def H6(info_rhd_file='', save_file=''):

    n = 64
    positions = np.zeros((n, 2))
    shank_ids = np.zeros((n,1),dtype=str)
    # shank_ids = []
    for i in range(n): # Two  shank probe
        x = (i % 2)*22.5 + (i//32)*250 # (um)
        y = (i % 32)*12.5
        positions[i] = x, y
        shank_ids[i]=str(i//32)
        # shank_ids.append((str(i//32)))

    # Associate native channel names (A-018 is the contact on the tip of the first shank, etc)
    mapping_to_device = np.array([
        18,50,16,26,17,14,54,32,44,30,12,52,27,28,11,9,46,22,24,20,15,3,48,43,13,
        5,45,23,19,47,1,21,39,55,37,4,61,51,57,33,62,38,35,63,56,0,25,53,7,42,36,
        41,6,8,58,60,34,2,31,40,59,29,49,10])

    if info_rhd_file:
    
        # To apply the custom channel map
        from src import rhdutilities 

        # Load the native channel names actually saved in this recording
        rhd_info = rhdutilities.load_file(info_rhd_file)
        rhd_info = rhd_info[0]
        amplifier_channels = rhd_info['amplifier_channels']
        n = len(amplifier_channels)
        native_channel_number = np.zeros((n,1))
        for i in range(n):
            native_channel_name = amplifier_channels[i]['native_channel_name']
            native_channel_number[i] = int(native_channel_name[2:])

        # Remove contacts not saved in this file
        mask = np.isin(mapping_to_device, native_channel_number)
        positions = positions[mask]      
        shank_ids = shank_ids[mask]
        mapping_to_device = mapping_to_device[mask] 

        # Re-number mapping to device to match .dat file indices
        mapping_to_device = np.array(np.argsort(np.argsort(mapping_to_device)))

    # Convert shank_id into an array of chars (there's probably a better way of doing this)
    shank_ids_final = []
    for i in range(n):
        shank_ids_final.append(shank_ids[i][0])
    shank_ids = shank_ids_final

    # Store positions and .dat file indices in the probe format
    probe = pi.Probe(ndim=2, si_units='um')
    probe.set_contacts(positions=positions, shapes='rect', 
        shape_params={'width': 11,'height':15}, shank_ids=shank_ids)
  
    # probe.create_auto_shape(probe_type='tip')
    probe.set_device_channel_indices(mapping_to_device)

    # Save to .json file on disk
    if save_file:
        pi.io.write_probeinterface(save_file, probe)

    return probe

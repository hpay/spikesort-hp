# spikesort-hp
My wrapper code for spike sorting with spike interface and/or kilosort

You should clone the github repositories:

- spikeinterface

- kilosort2.0

- ironclust


and place them in a flat directory structure, e.g.:

 ~\spikesort

    |--spikesort-hp
    
    |--kilosort2.0
    
    |--spikeinterface

Then navigate to the spikesort-hp\scripts folder and: 
- run_batch_kilosort.m to run kilo2.0 on a bunch of files in matlab using SC's custom settings
- runme.ipynb to run spikeinterface/whatever I'm working on currently 

Recommend to first run the built-in kilosort gui (kilosort.m in in the kilosort2.0 repo) to check kilosort installation

Then can run batch KS using our custom settings with run_batch_kilosort.m

I'm currently working on using spikeinterface to compare kilo2 with other spike sorters
Run the .ipynb or check_install_k2_iron.py to test install of spikeinterface and ability to call those sorters 

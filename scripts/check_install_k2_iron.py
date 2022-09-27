from pathlib import Path
import platform
import os
import shutil
import time

def _test():
    import spikeinterface.sorters as ss
    ss.installed_sorters()


def _create_recording():
    import spikeinterface.full as si
    rec, sorting = si.toy_example(num_segments=1, duration=100, seed=1, num_channels=16, num_columns=2)
    rec.save(folder='./results/toy_example_recording')

def _run_one_sorter_and_extract_wf(sorter_name):
    import spikeinterface.full as si
    rec = si.load_extractor('./results/toy_example_recording')
    sorting = si.run_sorter(sorter_name, rec, output_folder=f'results/{sorter_name}_output', verbose=True)
    si.extract_waveforms(rec, sorting, f'results/{sorter_name}_waveforms',
                         n_jobs=1, total_memory="10M", max_spikes_per_unit=500, return_scaled=False)


def run_tridesclous():
    _run_one_sorter_and_extract_wf('tridesclous')


def run_kilo2():
    _run_one_sorter_and_extract_wf('kilosort2-compiled')


def export_to_phy():
    import spikeinterface.full as si
    we = si.WaveformExtractor.load_from_folder("results/kilo2_waveforms")
    phy_folder = "./results/phy_example"
    si.export_to_phy(we, output_folder=phy_folder, verbose=False)


def open_phy():
    os.system("phy template-gui ./results/phy_example/params.py")

def _clean():
    # clean
    folders = [
        "kilo2_output", "kilo2_waveforms",
        "tridesclous_output", "tridesclous_waveforms",
        "spykingcircus_output", "spykingcircus_waveforms",
        "phy_example", "toy_example_recording"
    ]
    for folder in folders:
        if Path('results', folder).exists():
            shutil.rmtree(Path('results',folder))

if __name__ == '__main__':

    _test()

    _clean()
    _create_recording()

    steps = [
        ('Run kilo2', run_kilo2),
        ('Run tridesclous', run_tridesclous),
        ('Export to phy', export_to_phy),
        ('Open phy', open_phy),
    ]

    for label, func in steps:
        func()
        done = '...OK'
        print(label, done)


    _clean()

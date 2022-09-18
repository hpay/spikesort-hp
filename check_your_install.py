from pathlib import Path
import platform
import os
import shutil
import time

# If the "pip install -e ." in the spikeinterface folder doesn't work:
# import sys
# sys.path.insert(0, '../spikeinterface')


def check_import_si():
    import spikeinterface as si


def check_import_si_full():
    import spikeinterface.full as si


def _create_recording():
    import spikeinterface.full as si
    rec, sorting = si.toy_example(num_segments=1, duration=100, seed=1, num_channels=16, num_columns=2)
    rec.save(folder='./results/toy_example_recording')

def _run_one_sorter_and_exctract_wf(sorter_name):
    import spikeinterface.full as si
    rec = si.load_extractor('./results/toy_example_recording')
    sorting = si.run_sorter(sorter_name, rec, output_folder=f'results/{sorter_name}_output', verbose=False)
    si.extract_waveforms(rec, sorting, f'results/{sorter_name}_waveforms',
                         n_jobs=1, total_memory="10M", max_spikes_per_unit=500, return_scaled=False)


def run_tridesclous():
    _run_one_sorter_and_exctract_wf('tridesclous')


def run_spykingcircus():
    _run_one_sorter_and_exctract_wf('spykingcircus')


def export_to_phy():
    import spikeinterface.full as si
    we = si.WaveformExtractor.load_from_folder("results/tridesclous_waveforms")
    phy_folder = "./results/phy_example"
    si.export_to_phy(we, output_folder=phy_folder, verbose=False)


def open_phy():
    os.system("phy template-gui ./results/phy_example/params.py")


def _clean():
    # clean
    folders = [
        "tridesclous_output", "tridesclous_waveforms",
        "spykingcircus_output", "spykingcircus_waveforms",
        "phy_example", "toy_example_recording"
    ]
    for folder in folders:
        if Path('results', folder).exists():
            shutil.rmtree(Path('results',folder))


if __name__ == '__main__':

    _clean()
    _create_recording()

    steps = [
       # ('Import spikeinterface', check_import_si),
       # ('Import spikeinterface.full', check_import_si_full),
        ('Run tridesclous', run_tridesclous),
       # ('Open spikeinterface-gui', open_sigui),
        ('Export to phy', export_to_phy),
        ('Open phy', open_phy),
    ]

    if platform.system() == "Windows":
        pass
        #steps.insert(3, ('Run spykingcircus', run_spykingcircus))

    for label, func in steps:
        try:
            func()
            done = '...OK'
        except:
            done = '...Fail'
        print(label, done)

    _clean()

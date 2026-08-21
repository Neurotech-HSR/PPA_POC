# PPA_POC_SN — fMRI-derived EEG Functional Scouts for Subject-Specific rTMS

MATLAB/Brainstorm pipeline to derive **subject-specific EEG functional scouts** from fMRI-informed atlases, and to extract EEG features that discriminate **Silent Naming** from a contrast condition (**self-paced Counting**). The resulting ROIs are intended as candidate targets for subject-specific rTMS.

---

## Requirements

| Component | Notes |
|---|---|
| MATLAB | R2023 or later |
| [Brainstorm](https://neuroimage.usc.edu/brainstorm/) | Must be on the MATLAB path (handled in `pipeline_SN.m`) |
| Sample dataset | `yHC_Test.zip` — see [Data](#data) |

> Tested on Windows. Paths in the examples use Windows separators.

---

## Repository structure

```
PPA_POC_SN/
├── Config_SN.json              # Configuration file (paths + process switches)
├── pipeline_SN.m               # Main pipeline
├── SN_feature_extraction.m     # Feature extraction
├── SN_GUI_features.m           # GUI for ROI/feature inspection
├── data/                       # Data folder (target of `mainDir`)
└── Brainstorm_related/
    └── zipSubject/             # Place yHC_Test.zip here (do NOT unzip)
```

**Do not modify the folder organization.** All scripts resolve paths relative to this structure.

---

## Setup

### 1. Get the code

Download the repository `.zip` and unzip it, keeping the folder organization unchanged.

### 2. Get the sample subject

Download `yHC_Test.zip` from **https://drive.google.com/drive/folders/1WSh1cB4vg7vZC5B93z4whyiqrLNxy4EY?usp=sharing** and place it in:

```
PPA_POC_SN/Brainstorm_related/zipSubject/
```

⚠️ **Do not unzip this file.** Brainstorm imports it directly as an archive.

### 3. Configure paths

Open `Config_SN.json` in the main folder and set the `mainDir` field to the absolute path of the `data` folder:

```json
{
  "mainDir": "previousPath\\PPA_POC_SN\\data"
}
```

### 4. Set up the Brainstorm protocol

1. Open `pipeline_SN.m` and run only the command that adds Brainstorm to the MATLAB path, then launch Brainstorm.
2. Create a new protocol named **`SN_POC`**.
3. Import the subject archive: `File → Load protocol → Import subject from zip`, and select `yHC_Test.zip`.
4. Rename the imported subject to **`yHC_Test_T0_SN`**.

The subject name is used programmatically by the pipeline — the exact string matters.

---

## Running the pipeline

Run the three scripts in order:

| Step | Script | Output |
|---|---|---|
| 1 | `pipeline_SN.m` | Computes the files on which feature extraction is performed |
| 2 | `SN_feature_extraction.m` | Extracts EEG features per ROI |
| 3 | `SN_GUI_features.m` | GUI showing which ROIs display features that discriminate Silent Naming from the contrast condition |

---

## Data

`yHC_Test.zip` contains **already preprocessed MRI data** of a sample subject.

**Atlases.** Several functional atlases are included: some derived from CAT12 segmentation, others from in-house atlases. The **speech (MNI-cat12)** atlas is derived from the naming atlas available on [neurosynth.org](https://neurosynth.org/).

**Head mask surface.** Among the surface files there is a `head mask` surface, obtained with a **RevoScan 3D head scanner** and imported into Brainstorm. Electrode localization was performed on this surface and coregistered with the T1 data; the coregistered montage is stored in `channel.m`.

**Disabled processes.** In `Config_SN.json`, the MRI and fMRI import processes are **disabled by design**, since they have already been computed for this subject. Only the EEG files are imported; the channel file is added afterwards, already coregistered.

---

## Citation / Contact

<Add reference, license, and contact information here.>

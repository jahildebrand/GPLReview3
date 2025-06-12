# GPLReview3
Overview

GPLreview3 is an interactive MATLAB-based tool for reviewing, labeling, and validating acoustic detections. It displays concatenated spectrograms of detected sound events, allows users to listen to them, and assign or edit detection labels (true/false). The interface supports visualization of ship and airplane detections alongside animal detections and provides tools for audio filtering, playback at different speeds, and quantitative feedback on labeling progress.

GPLreview3 is intended for use on data processed through the GPL (Generalized Power Law) detector but can be adapted to similar workflows.
Getting Started

Launch the program: Run GPLreview3.m to open the GUI (GPLreview3.fig).

Load required files:

    Open a Detections file using the T icon (upper right). This provides event times.

    Open a WAV/XWAV file using the W icon. Select the first file in the folder of WAV files used for the GPL detections.

Both files are required: Detections for event times, WAV for raw audio data to generate spectrograms.

Main display: The window shows a concatenated spectrogram of all detections.
Navigation

    Forward button (>) — Steps forward to the next set of detections. Labeled Start Here.

    Backward button (<) — Steps backward.

    Plot button — Replots the current window.

Display Controls

    Brightness slider — Adjust spectrogram brightness.

    Whiten Spectrogram — Display pre-whitened spectrogram.

    Detection Labels — Shows detection boundaries as vertical lines and a colored bar:

        Green = true

        Red = false
        Default state: all false.

    Reset All True — Sets all detections to true. Requires confirmation.

    Airplane/Ship — Displays ship/airplane detections at the top of the spectrogram:

        Green/cyan = detection present

        Red = not present
        These labels cannot be edited in GPLreview.

Audio Playback

    Play Audio button — Activates crosshairs to select an interval to play:

        Click left → start

        Click right → end
        Two red lines will mark the selection.

    Filter Audio — After selecting the time range, click the lower-left and upper-right corners of the spectrogram. The vertical range of the crosshairs sets the filter frequency band.

    5X Speedup — Plays sound 5X faster (higher pitch and shorter playback time).

Plot Settings

    Plot Length (s) — Sets cumulative display time (default 75 s).

    Start Det — Sets starting detection number (default 1).
    This does not auto-update as you navigate.

    Progress — Percentage of detections reviewed is shown at the bottom.

    Frequency & FFT parameters — Read from input data and should not be changed.

Editing Detection Labels

    All Yes — Set all detections in the current window to true and advance to the next window.

    All No — Set all detections in the current window to false and advance.

    Some Yes — Use crosshairs to label a subset of detections as true.

    Some No — Use crosshairs to label a subset of detections as false.

Metrics Displayed

Displayed above the upper-left corner of the spectrogram:

    Total Time (min) — Total detection time.

    % False — Percentage of detection time labeled false.

    Total Airplane/Ship Time (min) — Time with airplane/ship detections.

    % False with Airplane/Ship — Percentage of airplane/ship time labeled false.

Saving Results

    Edits are saved in the Detections file (Labels variable):

        1 = true

        0 = false

    The Detections file also contains:

        parms — Detector parameters.

        Times — Start and stop times in samples and as MATLAB datenums.

Notes

    Detections longer than 1 second are split into multiple detections due to the GPL detector's behavior (merging closely spaced calls unless there is a sufficient gap).


# jBeamformers
Java beamformers for Ultrasound prebeamformed data module, called using matlab

1.    Download the ZIP file and extract it.

2.    Download an example pre-beamformed file from https://drive.google.com/drive/folders/1ZlI_EH8RFQw24tDwTIZ49dD3uClP7Dn4?usp=sharing.

        Data Source: This data is from Uliana, J. H., Sampaio, D. R. T., Fernandes, G. S. P., Brassesco, M. S., Nogueira-Barbosa, M. H., Carneiro, A. A. O., & Pavan, T. Z. (2020).
        Multiangle long-axis lateral illumination photoacoustic imaging using linear array transducer. Sensors (Switzerland), 20(14), 1–19. https://doi.org/10.3390/s20144052.

        Data Description: The examples represent a B-mode image and a Photoacoustic (PA) image of a human index finger.

 3.   Open RUN_ME.m in MATLAB (the code was developed in R2021a).

 4.   Select example_Bmode.mat to beamform the B-mode data, or example_PA.mat for the PA data.

 5.   The script will display the resulting images for the following beamforming techniques:

        Delay and Sum (DAS)

        Filtered Delay Multiply and Sum (fDMAS)

        Short-Lag Spatial Coherence (SLSC)

        Delay Multiply and Sum (DMAS)

<img width="1018" height="624" alt="image" src="https://github.com/user-attachments/assets/0ec2b305-3cf1-41f6-a89c-8055c0d066c6" />

<img width="994" height="611" alt="image" src="https://github.com/user-attachments/assets/c87d25db-2bd1-4c26-a39a-74a8344938d3" />


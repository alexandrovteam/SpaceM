def getPixSize(MFI):
    """Reads the pixel size in um from the Nikon Ti E microscope (NIS elements software).

    Args:
        MFI (str): path to Main Folder Input.

    Returns:
        pix_size (float): pixel size in um.

    """
    txt_file = codecs.open(MFI + '/Microscopy/postMALDI/out.txt', 'r', 'utf-16')
    for row in txt_file:
        if row.startswith('Calibration'):
            pix_size = float(row.strip().split()[2].replace(',', '.'))
        else:
            pix_size = 0.73  # Pixel size used in all experiments
    return pix_size
def read_check(
    read, 
    direction: bool, 
    barcode_tag: str = "CB", 
    barcode_len: int = 16, 
    seq_len: int = 98, 
    ignore_chro: list = None
):
    """
    Checks if a read is usable for analysis.
    """
    if ignore_chro is None:
        ignore_chro = ["MT", "M", "chrM"] # Standard defaults

    try:
        # Check if barcode exists and has correct length
        cb = read.get_tag(barcode_tag)
        if len(cb) != barcode_len:
            return 0, 0, 0, 0, 0
    except KeyError:
        return 0, 0, 0, 0, 0
    
    read_chro = read.reference_name
    read_start = read.reference_start
    read_end = read.reference_end
    read_strand = read.is_reverse
    
    # Extract sample ID safely (fallback if '.' is missing)
    try:
        sample_id = read.query_name.split('.')[0]
    except IndexError:
        sample_id = "Sample"

    # 1. Check Strand (Must match the direction we are currently analyzing)
    if direction != read_strand:
        return 0, 0, 0, 0, 0
     
    # 2. Check Chromosome Ignore List
    if read_chro in ignore_chro:
        return 0, 0, 0, 0, 0

    # 3. Normalizing Read Length
    # If read is too long, skip. If too short, extend it to standard length.
    current_len = read_end - read_start
    if current_len > seq_len:
        return 0, 0, 0, 0, 0
    elif current_len < seq_len:
        read_end = read_start + seq_len
    
    # Combine Sample ID and Barcode to ensure uniqueness across merged files
    cb_final = f"{sample_id}_{cb}"
    
    return read_chro, read_start, read_end, read_strand, cb_final

if __name__ == "__main__":
    pass

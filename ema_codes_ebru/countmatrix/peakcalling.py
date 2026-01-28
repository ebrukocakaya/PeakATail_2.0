import pysam as ps
import bisect
import time
from ema.countmatrix.peak import Peak
from ema.countmatrix.read import read_check
from ema.countmatrix.paswrite import matrix_write, pas_write

def peak_calling(
        direction: bool, 
        bedfilepath: str,
        matrixpath: str, 
        bamfile_dir: str, 
        default_threshold: int = 10, 
        merge_len: int = 50
    ):
    """
    Core function to detect PolyA sites from BAM file.
    Args:
        direction (bool): True for Negative strand, False for Positive.
        bedfilepath (str): Path to output BED file.
        matrixpath (str): Path to output matrix file.
        bamfile_dir (str): Path to input BAM file.
        default_threshold (int): Minimum coverage to consider a peak.
        merge_len (int): Distance to merge nearby peaks.
    """
    print(f"Processing BAM: {bamfile_dir} | Strand: {'Negative' if direction else 'Positive'}")
    
    bamfile = ps.AlignmentFile(bamfile_dir, 'rb')
    matrix = open(matrixpath, "w")
    bedfile = open(bedfilepath, "w")
    
    data_array = [] 
    signal = False
    chro = "1"
    l_end, i_end = 0, 0
    timercount = 0
    start_time = time.time()
    
    # Initialize Peak object
    peak = Peak(peak_strand=direction)
    i = 0 

    for read in bamfile:
        # Progress timer (every 1M reads)
        timercount += 1
        if timercount % 1000000 == 0:
            print(f"Processed {timercount} reads... ({time.time()-start_time:.2f}s)")
        
        # Check read validity
        # Note: read_check handles the logic of extracting tags/positions
        chro1, start1, end1, strand, cb = read_check(read=read, direction=direction)

        if chro1 == 0:
            continue
        
        # Chromosome Switch Logic (When moving from chr1 -> chr2)
        if chro1 != chro:
            if signal:
                pas_1, pas_2 = peak.pasfind()
                if pas_1 != 0:
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)
            
            elif len(peak.peak_list) != 0:
                pas_1, pas_2 = peak.pasfind()
                if pas_1 != 0:
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)
            
            signal = False
            # Reset Peak object for new chromosome
            peak = Peak(peak_start=0, peak_strand=direction, peak_list=[], cb_dict={}, last_peak_end=0)
            i_end, l_end, data_array, i = 0, 0, [], 0

        # Peak Processing Logic
        if signal:
            # Check the height at the threshold cut-off
            if len(data_array) >= default_threshold:
                l_end = data_array[-default_threshold]
                peak.cb_counting(cb=cb)

        # Check if peak has ended (current read start > last valid end)
        if signal and start1 > l_end: 
            signal = False
            peak.last_peak_end = l_end
            peak.peak_start = 0
            
        # Add current read end to the "pileup" array
        bisect.insort(data_array, end1)

        # Remove reads that have already ended (maintain sliding window)
        if start1 > i_end:
            slice_loc = bisect.bisect_left(data_array, start1)
            if signal:
                peak.peak_add(data_array=data_array, slice_loc=slice_loc)
            
            data_array = data_array[slice_loc:]
            if data_array:
                i_end = data_array[0]
            else:
                i_end = 0

        height = len(data_array)

        # Check if a new peak is starting
        if signal == False and height >= default_threshold:
            signal = True
            i += 1 

            # Check if this is a distinct peak or a merge-able one
            if start1 - peak.last_peak_end > merge_len and i != 1:
                # Previous peak is done, finalize it
                pas_1, pas_2 = peak.pasfind()

                if pas_1 != 0: 
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)

                # Start new peak instance
                peak = Peak(peak_start=start1, peak_strand=direction, peak_list=[], cb_dict={}, last_peak_end=0)
            else:
                # Merge with previous or start fresh
                peak.peak_start = start1 

        chro = chro1 

    bamfile.close()
    matrix.close()
    bedfile.close()

if __name__ == "__main__":
    # This block is for testing only; real execution comes from main.py
    pass

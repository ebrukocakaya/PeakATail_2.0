class Peak:
    pasnumber = 0

    def __init__(self, peak_list=None, peak_start=0, last_peak_end=0, peak_strand=True, cb_dict=None):
        """
        :param peak_list: List of [read_end, height]
        :param peak_start: Start position of the peak
        :param peak_strand: True (Negative?) / False (Positive?)
        :param cb_dict: Dictionary of CellBarcodes in this peak
        """
        self.peak_list = peak_list if peak_list is not None else []
        self.peak_start = peak_start
        self.last_peak_end = last_peak_end
        self.peak_strand = peak_strand
        self.cb_dict = cb_dict if cb_dict is not None else {}

    def peak_add(self, data_array: list, slice_loc: int):
        """Updates peak_list with new data from the sliding window."""
        array_len = len(data_array)
        for i in range(slice_loc):
            # Append [position, height]
            self.peak_list.append([data_array[i], array_len - i])

    def cb_counting(self, cb: str):
        """Increments count for a cell barcode."""
        if cb in self.cb_dict:
            self.cb_dict[cb] += 1
        else:
            self.cb_dict[cb] = 1
    
    def pasfind(self):
        """
        Finds the PolyA Site (PAS) within the peak.
        Logic: Finds the max height, calculates 5% threshold, and finds the range.
        """
        try:
            if not self.peak_list:
                return 0, 0

            max_height = max(self.peak_list, key=lambda x: x[1])[1]
            
            # Threshold Check (5% of max height must be significant)
            if max_height <= 20:
                return 0, 0
            
            threshold = max_height * 0.05
            pas_1, pas_2, pas_cov = 0, 0, 0

            # Logic differs based on strand direction
            if not self.peak_strand:
                # Positive Strand Logic
                for item in reversed(self.peak_list[:-1]):
                    pos, cov = item[0], item[1]
                    if cov >= threshold and pas_1 == 0:
                        pas_1, pas_cov = pos, cov
                    elif pos <= pas_1 and cov > pas_cov:
                        pas_2 = pos
                        return pas_1, pas_2
            else:
                # Negative Strand Logic
                for item in self.peak_list:
                    pos, cov = item[0], item[1]
                    if cov >= threshold and pas_1 == 0:
                        pas_1, pas_cov = pos, cov
                    elif pos >= pas_1 and cov > pas_cov:
                        pas_2 = pos
                        return pas_1, pas_2
            
            return 0, 0
        except Exception as e:
            # print(f"Error in pasfind: {e}") # Debugging if needed
            return 0, 0

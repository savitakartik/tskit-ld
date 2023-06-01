import numpy as np
import tskit
import time

class LdInfo:
    def __init__(self, ts, chr):
        self.ts = ts
        self.chr = chr

    def get_single_mutation_sites(self):
        sites_with_one_mut = np.where(
            np.bincount(self.ts.mutations_site, minlength=self.ts.num_sites) == 1
        )[0]
        return sites_with_one_mut

    def compute_ld_with_time(self, sites_a, sites_b):
        ldcalc = tskit.LdCalculator(self.ts)
        r2 = []
        compute_times = []
        for a, b in zip(sites_a, sites_b):
            start = time.perf_counter()
            r2.append(ldcalc.r2(a, b))
            compute_times.append(time.perf_counter() - start)
        return r2, compute_times

    def compute_distance_between_sites(self, sites_a, sites_b):
        distances = []
        for a, b in zip(sites_a, sites_b):
            distances.append(abs(self.ts.sites_position[a] - self.ts.sites_position[b]))
        return distances

    def return_random_sites_in_range(self, sites_a, max_dist=200_000):
        sites_with_one_mut = self.get_single_mutation_sites()
        sites_pos_with_one_mut = np.take(self.ts.sites_position, sites_with_one_mut)
        sites_b = []
        for a in sites_a:
            a_pos = self.ts.sites_position[a]
            ld_range = np.where(
                (sites_pos_with_one_mut >= a_pos - max_dist)
                & (sites_pos_with_one_mut <= a_pos + max_dist)
            )[0]
            rand_ind = np.random.choice(ld_range)
            sites_b.append(sites_with_one_mut[rand_ind])
        return sites_b



'''
A place to store some prototype or experimental code, 
which may be garbage or useful in someway elsewhere.
'''



# removed from class epdsConstructor:
# -----------------------------------------------------------------------------
#
    def exhaustive_build(self, ii, list_empCpds):
        '''
        exhaustive_search=False, exclude_singletons=True
                # the now remaining are to be assigned as singletons
        if not exclude_singletons:
            remaining_peaks = [P for P in remaining_peaks if P['id_number'] not in found2]
            for P in remaining_peaks:
                ii += 1
                list_empCpds.append(
                    {'id': ii, 'list_peaks': [(P['id_number'], 'anchor')]}
                )

        '''
        pass

    def __extend_empCpds_by_adducts(self, seed_list_empCpds, list_peaks, adduct_patterns, mz_tolerance_ppm=5):
        '''
        Not used now.
        Search list_peaks for adducts that fit patterns relative to anchor ions in existing empCpds.
        Co-elution required.
        This is not done in initial empCpd construction to reduce errors in assigning peaks to empCpds.

        in progress -
        '''
        
        mztree = build_centurion_tree(list_peaks)
        for EPD in seed_list_empCpds:
            P1, relation = EPD['list_peaks'][0]
            for adduct in adduct_patterns:
                # (1.0078, 'H'), (21.9820, 'Na/H'), ...
                tmp = find_all_matches_centurion_indexed_list( P1['mz'] + adduct[0], mztree, mz_tolerance_ppm )
                for P2 in tmp:
                    if is_coeluted(P1, P2):
                        EPD['list_peaks'].append( (P2['id_number'], relation +','+ adduct[1]) )


        return seed_list_empCpds


    def __consolidate_epds(self, epd_peaks):
        '''
        Merge redundant empirical compounds (epds), which can arise when diff isotopes have same adducts,
        by 1) removing an epd if it's a subset of another epd;
           2) resolving each peak to only one epd.
        epd_peaks = [[(182, 'anchor'), (191, '13C/12C'), (205, '18O/16O')], ...]
        '''
        ordered_set_peak_ids = [set([x[0] for x in P]) for P in epd_peaks]
        def _choose(jj, ii, ordered_set_peak_ids=ordered_set_peak_ids, epd_peaks=epd_peaks):
            if ordered_set_peak_ids[ii].issubset(ordered_set_peak_ids[jj]):
                return jj
            elif ordered_set_peak_ids[jj].issubset(ordered_set_peak_ids[ii]):
                return ii
            else:           # warning for now; should choose the more common ion
                print("Unresolved peak between ", ii, jj)
                #
                # to update
                #
                return jj

        p2e = {}
        for ii in range(len(ordered_set_peak_ids)):
            for p in ordered_set_peak_ids[ii]:
                if p not in p2e:
                    p2e[p] = ii
                else:
                    p2e[p] = _choose(p2e[p], ii)

        valid_epds = set(p2e.values())
        # this does not completely resolve peak ownership
        return [epd_peaks[ii] for ii in valid_epds]



def filter_overlapping_orfs(orfs):
    """
    Remove overlapping ORFs and keep the longest one.
    """
    filtered_orfs = []
    orfs = sorted(orfs, key=lambda x: (x[0], x[1]))
    for start, end, orf_seq in orfs:
        if not filtered_orfs or start > filtered_orfs[-1][1]:
            filtered_orfs.append((start, end, orf_seq))
        else:
            if (end - start) > (filtered_orfs[-1][1] - filtered_orfs[-1][0]):
                filtered_orfs[-1] = (start, end, orf_seq)
    return filtered_orfs

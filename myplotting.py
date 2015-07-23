import numpy as np


def delete_overlapping_xlabels(fig, ax):

    ''' Deletes overlapping xtick labels'''

    import numpy as np

    fig.canvas.draw()
    major_labels = ax.get_xticklabels(minor=False)
    new_major_text = [item.get_text() for item in major_labels]

    print new_major_text

    old_box = None
    for i, label in enumerate(major_labels):
        new_box = label.get_window_extent().get_points()[:, 0]

        if old_box is not None:
            if ((new_box[0] < old_box[1]) & \
                (new_box[0] > old_box[0])):
                new_major_text[i] = ''
                new_box = [np.nan, np.nan]

        old_box = new_box

    ax.set_xticklabels(new_major_text)

    # minor ticks
    fig.canvas.draw()
    minor_labels = ax.get_xticklabels(minor=True)
    new_minor_text = [item.get_text() for item in minor_labels]

    print new_minor_text

    old_box = None
    for i, label in enumerate(minor_labels):
        new_box = label.get_window_extent().get_points()[:, 0]

        print new_box, new_minor_text[i]

        if old_box is not None:
            if (((new_box[0] > old_box[0]) & \
                 (new_box[0] < old_box[1])) | \
                ((new_box[1] > old_box[0]) & \
                 (new_box[1] < old_box[1]))):
                print 'bounded'
                new_minor_text[i] = ''
                new_box = [np.nan, np.nan]

        old_box = new_box

    ax.set_xticklabels(new_minor_text, minor=True)

    # minor ticks against major
    fig.canvas.draw()
    minor_labels = ax.get_xticklabels(minor=True)
    major_labels = ax.get_xticklabels(minor=False)
    new_minor_text = [item.get_text() for item in minor_labels]

    for i, minor_label in enumerate(minor_labels):
        box_minor = minor_label.get_window_extent().get_points()[:, 0]

        for major_label in major_labels:
            box_major = \
                    major_label.get_window_extent().get_points()[:, 0]

            if (((box_minor[0] > box_major[0]) & \
                 (box_minor[0] < box_major[1])) | \
                ((box_minor[1] > box_major[0]) & \
                 (box_minor[1] < box_major[1]))):
                new_minor_text[i] = ''

    ax.set_xticklabels(new_minor_text, minor=True)

    return ax


def delete_overlapping_xlabels2(fig, ax):

    ''' Deletes overlapping xtick labels'''

    import numpy as np

    fig.canvas.draw()
    major_labels = ax.get_xticklabels(minor=False)
    new_major_text = [item.get_text() for item in major_labels]

    print new_major_text

    old_box = None
    for i in xrange(len(major_labels)):
        label = major_labels[i]
        new_box = label.get_window_extent().get_points()[:, 0]

        if old_box is not None:
            for j in xrange(i):
                old_box = major_labels[j].get_window_extent().get_points()[:, 0]
                if new_box[1] > old_box[0]:
                    new_major_text[i] = ''
                    new_box = [np.nan, np.nan]

        old_box = new_box

    print new_major_text
    ax.set_xticklabels(new_major_text)

    # minor ticks
    fig.canvas.draw()
    minor_labels = ax.get_xticklabels(minor=True)
    new_minor_text = [item.get_text() for item in minor_labels]

    print new_minor_text

    old_box = None
    for i in xrange(len(minor_labels)):
        label = minor_labels[i]
        new_box = label.get_window_extent().get_points()[:, 0]

        if old_box is not None:
            # check
            for j in xrange(i):
                old_box = minor_labels[j].get_window_extent().get_points()[:, 0]
                if new_box[1] > old_box[0]:
                    new_minor_text[i] = ''
                    new_box = [np.nan, np.nan]

        old_box = new_box

    ax.set_xticklabels(new_minor_text, minor=True)

    # minor ticks against major
    fig.canvas.draw()
    minor_labels = ax.get_xticklabels(minor=True)
    major_labels = ax.get_xticklabels(minor=False)
    new_minor_text = [item.get_text() for item in minor_labels]

    for i, minor_label in enumerate(minor_labels):
        box_minor = minor_label.get_window_extent().get_points()[:, 0]

        for major_label in major_labels:
            box_major = \
                    major_label.get_window_extent().get_points()[:, 0]

            if (((box_minor[0] > box_major[0]) & \
                 (box_minor[0] < box_major[1])) | \
                ((box_minor[1] > box_major[0]) & \
                 (box_minor[1] < box_major[1]))):
                new_minor_text[i] = ''

    ax.set_xticklabels(new_minor_text, minor=True)

    return ax


def delete_overlapping_xlabels3(fig, ax):

    ''' Deletes overlapping xtick labels'''

    import numpy as np

    fig.canvas.draw()
    major_labels = ax.get_xticklabels(minor=False)
    new_major_text = [item.get_text() for item in major_labels]
    bboxes = [label.get_window_extent() for label in major_labels]

    bbox_overlaps = check_overlaps(bboxes)
    overlaps = any(bbox_overlaps)

    while overlaps:
        i = np.argmax(bbox_overlaps)
        new_major_text[i] = ''
        bboxes[i].set_points(np.array([[np.nan, np.nan], [np.nan, np.nan]]))
        bbox_overlaps = check_overlaps(bboxes)
        overlaps = any(bbox_overlaps)

    ax.set_xticklabels(new_major_text)


def check_overlaps(bboxes):

    '''
    takes list of bboxes
    returns a list of how many times each bbox overlaps with other bboxes
    '''

    overlaps = [0] * len(bboxes)
    for i, box in enumerate(bboxes):
        for other_box in bboxes:
            if (box != other_box):
                #overlaps[i] += box.overlaps(other_box)
                overlaps[i] += check_overlaps(box, other_box)
    return overlaps


def check_overlaps(bbox1, bbox2):
    """
    Returns True if this bounding box overlaps with the given
    bounding box *other*. OR IF NANA!
    """
    ax1, ay1, ax2, ay2 = bbox1._get_extents()
    bx1, by1, bx2, by2 = bbox2._get_extents()

    if ax2 < ax1:
        ax2, ax1 = ax1, ax2
    if ay2 < ay1:
        ay2, ay1 = ay1, ay2
    if bx2 < bx1:
        bx2, bx1 = bx1, bx2
    if by2 < by1:
        by2, by1 = by1, by2

    if np.isnan(np.sum(extents1)) or np.isnan(np.sum(extents2)):
        return False

    return not ((bx2 < ax1) or
                (by2 < ay1) or
                (bx1 > ax2) or
                (by1 > ay2))

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):

    ''' Truncates a matplolib.colors colormap to a smaller range.

    Parameters
    ----------
    cmap : matplotlib.pyplot.cm
        Colormap
    minval : float
        Lower value to truncate.
    maxval : float
        Upper value to truncate
    n : int
        Number of discrete samples of colormap between minval and maxval.

    Returns
    -------
    new_cmap : matplotlib.pyplot.cm
        Truncated colormap

    '''

    import matplotlib.colors as colors

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))

    return new_cmap




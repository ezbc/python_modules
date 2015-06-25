import numpy as np
import matplotlib.pyplot as plt
import myplotting as myplt
from matplotlib.ticker import ScalarFormatter

def test_textoverlap_simple():
    plt.close(); plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(range(10), range(10))
    testing = True
    if testing:
    # Mess things up
        fig.canvas.draw()
        labels = [item.get_text().zfill(10) for item in ax.get_xticklabels()]
        ax.set_xticklabels(labels)

    fig.canvas.draw()

    ax = myplt.delete_overlapping_xlabels3(fig,ax)

    #plt.show()

def test_textoverlap_minorticks():

    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    minorLocator   = MultipleLocator(5)



    plt.close(); plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(range(100), range(100))
    testing = True
    # Mess things up
    if testing:
        fig.canvas.draw()
        labels = [item.get_text().zfill(3) for item in \
                ax.get_xticklabels(minor=True)]
        ax.set_xticklabels(labels, minor=True)
        fig.canvas.draw()
        #ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_minor_formatter(FormatStrFormatter('000000'))
        #plt.setp(ax.get_xticklabels(minor=True), visible=True)

    fig.canvas.draw()

    #ax = myplt.delete_overlapping_xlabels3(fig,ax)

    plt.show()

def test_bbox():

    import numpy as np
    import matplotlib.pyplot as plt


    plt.close(); plt.clf()
    fig, ax = plt.subplots(1)
    ax.plot(range(10), range(10))
    fig.canvas.draw()
    bbox1 = ax.get_xticklabels()[0].get_window_extent()
    bbox2 = ax.get_xticklabels()[1].get_window_extent()
    print('Overlap before nan values', bbox1.overlaps(bbox2))
    overlap = bbox2.set_points(np.array([[np.nan, np.nan], [np.nan, np.nan]]))
    print('Overlap after nan values', bbox1.overlaps(bbox2))
    assert not overlap

    '''
        import numpy as np
        import matplotlib.pyplot as plt


        plt.close(); plt.clf()
        fig, ax = plt.subplots(1)
        ax.plot(range(10), range(10))
        fig.canvas.draw()
        bbox1 = ax.get_xticklabels()[0].get_window_extent()
        bbox2 = ax.get_xticklabels()[1].get_window_extent()
        print('Overlap before nan values', bbox1.overlaps(bbox2))
        overlap = bbox2.set_points(np.array([[np.nan, np.nan], [np.nan, np.nan]]))
        print('Overlap after nan values', bbox1.overlaps(bbox2))
    '''


"""
tree_queries.py

Miscellaneous functions that extract data from an MDSplus tree.

In this initial version, the functions contained within this file are
specifically tailored to the implementation of the real-time system on the
LHD experiment. As such, they may not work for other implementations in which
the tree structure may not be exactly the same.
"""

import MDSplus as mds
import numpy as np

TREE_NAME='rt_mpts'
DIG_TREE_NAME='gtchilin'
N_SAMPLES=606
IND_DIG_CARD  = 2
IND_DIG_INPUT = 4

def getVoltageData(shot, shelf, poly, nSamples=N_SAMPLES, \
                   treeName=TREE_NAME, digTreeName=DIG_TREE_NAME):
    """
    Obtains raw polychromator signals in volts for a given shot, spatial
    volume, and polychromator channel.

    Parameters
    ----------
        shot: integer
            Shot number
        shelf: integer
            ID number of the "shelf" (spatial volume) in the tree
        poly: integer
            ID number of the polychromator channel in the subtree for 
            the spatial volume denoted by `shelf`
        nSamples: integer (optional)
            Number of samples collected (time points) for each laser pulse
        treeName: string (optional)
            Name associated with the tree files for experimental data. Note 
            that the location of the tree files must also be stored in an 
            environment variable with the name `[treeName]_path`.
        digTreeName: string (optional)
            Name associated with the tree files for digitizer calibration data. 
            Note that the location of the tree files must also be stored in an 
            environment variable with the name `[digTreeName]_path`.

    Returns
    -------
        signal: 2D double array
            Array containing the requested signal in volts. Each row corresponds
            to a successive laser pulse within the shot. Each column corresponds
            to a successive time point in its respective row (laser pulse).
    """

    shelf_str = 'shelf_%02d' % (shelf)
    poly_str  = 'poly_%d' % (poly)

    # Extract the raw signal from the signal/shot tree
    trSig = mds.Tree(treeName, shot, 'ReadOnly')
    dataNode = trSig.getNode('data').getNode(shelf_str).getNode(poly_str)
    rawSig = np.reshape(dataNode.getData().data(), (-1,nSamples)).T
    
    # Determine which digitizer card/input the signal was obtained from
    cardStr = trSig.getNode('mapping').getNode('shelf').getNode(shelf_str) \
                   .getNode('poly').getNode(poly_str).getNode('daq') \
                   .getData().data()   
    cardSubstrs = cardStr.split('_')
    cardNum = int(cardSubstrs[IND_DIG_CARD]) + 1
    cardInp = int(cardSubstrs[IND_DIG_INPUT]) + 1

    # Retreive calibration factors from the digitizer tree
    trDig = mds.Tree(digTreeName, -1, 'ReadOnly')
    inpNode = trDig.getNode('TS_S3316_%02d' % (cardNum)) \
                   .getNode('INPUT_%02d' % (cardInp))
    gain   = inpNode.getNode('gain_span').getData().data()
    offset = inpNode.getNode('offset_corr').getData().data()

    return ((rawSig - offset)/gain - 1.)*10.

def getRawData(shot, shelf, poly, nSamples=N_SAMPLES, treeName=TREE_NAME):
    """
    Obtains raw polychromator signals in terms of the unscaled digitizer
    signals for a given shot, spatial volume, and polychromator channel.
    This query is independent of the digitizer calibration data.

    Parameters
    ----------
        shot: integer
            Shot number
        shelf: integer
            ID number of the "shelf" (spatial volume) in the tree
        poly: integer
            ID number of the polychromator channel in the subtree for 
            the spatial volume denoted by `shelf`
        nSamples: integer (optional)
            Number of samples collected (time points) for each laser pulse
        treeName: string (optional)
            Name associated with the tree files for experimental data. Note 
            that the location of the tree files must also be stored in an 
            environment variable with the name `[treeName]_path`.

    Returns
    -------
        signal: 2D double array
            Array containing the requested signal in terms of digitizer values,
            all of which are positive integers. Each row corresponds to a 
            successive laser pulse within the shot. Each column corresponds to
            a successive time point in its respective row (laser pulse).

    """

    tr = mds.Tree(treeName, shot, 'ReadOnly')
    shelf_str = 'shelf_%02d' % (shelf)
    poly_str  = 'poly_%d' % (poly)
    dataNode = tr.getNode('data').getNode(shelf_str).getNode(poly_str)
    return  np.reshape(dataNode.getData().data(), (-1,nSamples)).T



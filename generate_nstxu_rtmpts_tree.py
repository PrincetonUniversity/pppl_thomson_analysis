"""
generate_nstxu_rtmpts_tree.py

Contains functions for generating an MDSplus model tree file for storing
calibration and experimental data from the NSTX-U real-time multi-point 
Thomson scattering (rt-mpts) system.

Usage: python generate_nstxu_rtmpts_tree.py [tree_name]

Note: to create a tree with the name [tree_name], it is first necessary to 
    create an environment variable [tree_name]_path containing the location 
    where the tree files are to be stored.

Maintained by: K.C. Hammond
Last update:   2021-06-17
Contact:       khammond@pppl.gov
"""

import MDSplus as mds
import numpy as np
import sys
import os

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

# Quantities
n_los              =  42   # number of lines of sight
n_polychromators   =   6   # number of polychromator outputs per line of sight
n_realTimeChannels =   8   # number of real-time channels
n_lasers           =   2   # number of lasers in use


# Base names for tree nodes
losName                    = 'shelf'
polychromatorName          = 'poly'
laserName                  = 'laser'
realTimeChannelName        = 'chan'
losShortName               = 'sh'
polychromatorShortName     = 'po'
realTimeChannelShortName   = 'ch'

# Root notes in the tree
rootNodes = ['calibration', 'daq', 'data', 'logs', 'mapping', 'rt', \
             'trigtime', 'eventtimes']

# Names for properties associated with each laser line-of-sight
#     radius  = radial location of the diagnosed volume
#     dRadius = uncertainty in radius
#     theta   = scattering angle (degrees)
#     cg      = coefficients for a polynomial expansion used for T_e guess...
#     wv_all  = wavelengths (angstroms) for polychromator transmissivity spectra
#     fb      = wavelength-dependent calib. factor for input to polychromators
#     dfb     = error in fb
losPropNames = ['radius', 'dRadius', 'theta', 'cg', 'wv_all', 'fb', 'dfb']

# Names of LOS and laser properties (nf = density factor; dnf = uncertainty)
losLaserPropNames = ['nf', 'dnf']

# Mapping names for lasers
signalMappingPropNames =  ['active', 'rt_active', 'rt_output']

laserSignalNames = ['R_st', 'R_sp']

laserMappingPropNames = ['active', 'rt_active', 'daq']

# Names for properties associated with polychromator output signals
#     qt_lambda: quantum transmissivity as a function of wavelength
#     dqt_lambda: uncertainty in qt_lambda
polychromatorPropNames = ['M', 'qt_lambda', 'dqt_lambda', 'cfb', 'dcfb', \
                          'cf1', 'cf2', 'gfast', 'gslow']

polychromatorMappingPropNames = ['active', 'rt_active', 'daq']

# Mapping names for real-time channels
channelMappingPropNames = ['active', 'rt_active', 'rt_output']

#-------------------------------------------------------------------------------
# Main tree-building function
#-------------------------------------------------------------------------------
def make_tree(treeName):

    # Verify that a suitable environment variable exists for the tree
    if not (treeName + '_path') in os.environ:
        raise Exception('Environment variable ' + treeName + '_path ' \
                        'for tree ' + treeName + ' not found')

    tree = mds.Tree(treeName, -1, 'new')
    topNode = tree.getDefault()

    build_node_calibration(topNode)
    build_node_mapping(topNode)
    build_node_daq(topNode)
    build_node_data(topNode)
    build_node_logs(topNode)
    build_node_rt(topNode)

    tree.write()


#-------------------------------------------------------------------------------
# Calibration subtree
#-------------------------------------------------------------------------------
def build_node_calibration(topNode):
    """
    Stores calibration data for the lines of sight, laser energy signals,
    and polychromator signals. 
    
    Currently does not contain calibration data relevant for the data 
    acquisition cards.
    """

    calibNode = topNode.addNode('calibration', 'structure')
    calibNode.setCompressOnPut(True)

    # List of source files for calibration data
    fileNode = calibNode.addNode('input_files', 'text')
    fileNode.setWriteOnce(False)
    fileNode.setCompressOnPut(True)

    # Calibration data for lines of sight
    losNode = calibNode.addNode(losName, 'structure')

    for i in range(n_los):

        numStr = '_%%0%dd' % np.ceil(np.log10(n_los))
        losNode_i = losNode.addNode(losName + numStr % (i+1), 'structure')

        # Geometric properties of the lines of sight
        for losPropName in losPropNames:
            losPropNode = losNode_i.addNode(losPropName, 'numeric')
            losPropNode.setWriteOnce(False)
            losPropNode.setCompressOnPut(True)

        # Laser-specific properties
        for losLaserPropName in losLaserPropNames:
            for j in range(n_lasers):
                name_j = laserName + '_%d_' % (j+1) + losLaserPropName
                losLaserPropNode = losNode_i.addNode(name_j, 'numeric')
                losLaserPropNode.setWriteOnce(False)
                losLaserPropNode.setCompressOnPut(True)

        # Polychromator properties
        polyNode = losNode_i.addNode(polychromatorName, 'structure')
        for j in range(n_polychromators):
            polyNode_j = polyNode.addNode( \
                             polychromatorName + '_%d' % (j+1), 'structure')
            for polychromatorPropName in polychromatorPropNames:
                polyPropNode = \
                    polyNode_j.addNode(polychromatorPropName, 'numeric')
                polyPropNode.setWriteOnce(False)
                polyPropNode.setCompressOnPut(True)

            # there may be further properties for POLYCHROMATORSIGNALSLASER

    # Calibration data for lasers
    laserNode = calibNode.addNode(laserName, 'structure')

    for i in range(n_lasers):
        laserName_i = laserName + '_%d' % (i+1)
        laserNode_i = laserNode.addNode(laserName_i, 'structure')

        # Nodes for properties inherent to the laser
        for laserSignalName in laserSignalNames:
            laserSignalNode = laserNode_i.addNode(laserSignalName, 'numeric')
            laserSignalNode.setWriteOnce(False)
            laserSignalNode.setCompressOnPut(True)

        # Cross-references to properties specific to the lines of sight
        for losLaserPropName in losLaserPropNames:
            losPropNode = laserNode_i.addNode(losLaserPropName, 'structure')
            for j in range(n_los):
                numStr = '_%%0%dd' % np.ceil(np.log10(n_los))
                losName_j = losName + numStr % (j+1)
                losNode_j = losPropNode.addNode(losName_j, 'numeric')

                # Cross-reference to corresponding node in losNode subtree
                losNode_j.putData( \
                    losNode.getNode(losName_j).getNode( \
                        laserName_i + '_' + losLaserPropName))

        # Note: code exists to add structures with names in 
        # "POLYCHROMATORSIGNALSLASER", but this is an empty array

    # Shutter configuration (NSTX-specific)
    shutter_conf_node = laserNode.addNode('shutter_conf', 'text');
    shutter_conf_node.setWriteOnce(False)
    shutter_conf_node.setCompressOnPut(True)

#-------------------------------------------------------------------------------
# Mapping subtree
#-------------------------------------------------------------------------------
def build_node_mapping(topNode):

    mappingNode = topNode.addNode('mapping', 'structure')
    mappingNode.setCompressOnPut(True)

    laserNode = mappingNode.addNode(laserName, 'structure')

    # Mapping to the analog inputs for each of the laser ID signals
    for i in range(n_lasers):
        laserName_i = laserName + '_%d' % (i+1)
        laserNode_i = laserNode.addNode(laserName_i, 'structure')
        for propName in laserMappingPropNames:
            propNode = laserNode_i.addNode(propName, 'text')
            propNode.setWriteOnce(False)
            propNode.setCompressOnPut(True)

    # Mapping to the analog inputs for the laser energy signal
    laserNode_energy = laserNode.addNode(laserName + '_energy', 'structure')
    for propName in laserMappingPropNames:
        propNode = laserNode_energy.addNode(propName, 'text')
        propNode.setWriteOnce(False)
        propNode.setCompressOnPut(True)

    losNode = mappingNode.addNode(losName, 'structure')

    # Mapping information for data from each laser line-of-sight
    for i in range(n_los):
        numStr = '_%%0%dd' % np.ceil(np.log10(n_los))
        losNode_i = losNode.addNode(losName + numStr % (i+1), 'structure')
        for propName in signalMappingPropNames:
            propNode = losNode_i.addNode(propName, 'text')
            propNode.setWriteOnce(False)
            propNode.setCompressOnPut(True)

        # Mapping to the analog inputs for each polychromator output
        polyNode = losNode_i.addNode(polychromatorName, 'structure')
        for j in range(n_polychromators):
            polyNode_j = polyNode.addNode( \
                             polychromatorName + '_%d' % (j+1), 'structure')
            for propName in polychromatorMappingPropNames:
                propNode = polyNode_j.addNode(propName, 'text')
                propNode.setWriteOnce(False)
                propNode.setCompressOnPut(True)

    chanNode = mappingNode.addNode(realTimeChannelName, 'structure')
    
    # Mapping information for real-time channels
    #     must be populated with cross-references to corresponding nodes 
    #     under losNode during the general model tree populating procedures
    for i in range(n_realTimeChannels):

        # To contain cross-references to info on the los assigned to the channel
        numStr = '_%%0%dd' % np.ceil(np.log10(n_realTimeChannels))
        chanNode_i = chanNode.addNode( \
                         realTimeChannelName + numStr % (i+1), 'structure')
        for propName in signalMappingPropNames:
            propNode = chanNode_i.addNode(propName, 'text')
            propNode.setWriteOnce(False)
            propNode.setCompressOnPut(True)
        
        # To contain cross-references to polychromator assigned to the channel
        polyNode = chanNode_i.addNode(polychromatorName, 'structure')
        for j in range(n_polychromators):
            polyNode_j = polyNode.addNode( \
                             polychromatorName + '_%d' % (j+1), 'structure')
            for propName in polychromatorMappingPropNames:
                propNode = polyNode_j.addNode(propName, 'text')
                propNode.setWriteOnce(False)
                propNode.setCompressOnPut(True)


    # TODO: add Struck card mapping information


#-------------------------------------------------------------------------------
# DAQ subtree
#
# Currently empty.
#-------------------------------------------------------------------------------
def build_node_daq(topNode):

    daqNode = topNode.addNode('daq', 'structure')
    daqNode.setCompressOnPut(True)

    # So far not filled

#-------------------------------------------------------------------------------
# Data subtree
# 
# Raw data from the DAQ cards as well as outputs from the real-time processing
# code will be stored here.
#-------------------------------------------------------------------------------
def build_node_data(topNode):

    dataNode = topNode.addNode('data', 'structure')
    dataNode.setCompressOnPut(True)

    for i in range(n_los):
        numStr = '_%%0%dd' % np.ceil(np.log10(n_los))
        losNode = dataNode.addNode(losName + numStr % (i+1), 'structure')
        for j in range(n_polychromators):
            sigNode = losNode.addNode( \
                          polychromatorName + '_%d' % (j+1), 'signal')
            sigNode.setWriteOnce(False)
            sigNode.setCompressOnPut(True)

    for i in range(n_lasers):
        sigNode = dataNode.addNode(laserName + '_%d' % (i+1), 'signal')
        sigNode.setWriteOnce(False)
        sigNode.setCompressOnPut(True)

    sigNode = dataNode.addNode(laserName + '_energy', 'signal')
    sigNode.setWriteOnce(False)
    sigNode.setCompressOnPut(True)

#-------------------------------------------------------------------------------
# Logs subtree
# 
# Currently empty.
#-------------------------------------------------------------------------------
def build_node_logs(topNode):

    logsNode = topNode.addNode('logs', 'structure')
    logsNode.setCompressOnPut(True)

    # So far not filled

#-------------------------------------------------------------------------------
# RT subtree
#
# Appears to contain much of the information already available in the 
# Mapping subtree... could this be redundant?
#-------------------------------------------------------------------------------
def build_node_rt(topNode):

    rtNode = topNode.addNode('rt', 'structure')
    rtNode.setCompressOnPut(True)

    # Subnode for laser information
    laserNode = rtNode.addNode(laserName, 'structure')
    for i in range(n_lasers):

        laserNode_i = laserNode.addNode(laserName + '_%d' % (i+1), 'structure')

        for laserSignalName in laserSignalNames:
            sigNode = laserNode_i.addNode(laserSignalName, 'numeric')
            sigNode.setWriteOnce(False)
            sigNode.setCompressOnPut(True)

        for propName in laserMappingPropNames:
            propNode = laserNode_i.addNode(propName, 'text')
            propNode.setWriteOnce(False)
            propNode.setCompressOnPut(True)

        for losPropName in losLaserPropNames:
            losNode = laserNode_i.addNode(losPropName, 'structure')
            for j in range(n_realTimeChannels):
                numStr = '_%%0%dd' % np.ceil(np.log10(n_realTimeChannels))
                chanNode = losNode.addNode( \
                               realTimeChannelName + numStr % (j+1), 'numeric')
                chanNode.setWriteOnce(False)
                chanNode.setCompressOnPut(True)

        # Florian's script (line 487) also adds a "Poly" substructure here, but 
        # this doesn't appear in the tree on nstx server (or the rt_mpts tree)

    laserEnergyNode = laserNode.addNode(laserName + '_energy', 'structure')
    for propName in laserMappingPropNames:
        propNode = laserEnergyNode.addNode(propName, 'text')
        propNode.setWriteOnce(False)
        propNode.setCompressOnPut(True)

    shutterNode = laserNode.addNode('shutter_conf', 'text')
    shutterNode.setWriteOnce(False)
    shutterNode.setCompressOnPut(True)

    # Subnode for real-time channel information
    chanNode = rtNode.addNode(realTimeChannelName, 'structure')
    for i in range(n_realTimeChannels):

        numStr = '_%%0%dd' % np.ceil(np.log10(n_realTimeChannels))
        chanNode_i = chanNode.addNode( \
                         realTimeChannelName + numStr % (i+1), 'structure')

        for propName in losPropNames:
            propNode = chanNode_i.addNode(propName, 'numeric')
            propNode.setWriteOnce(False)
            propNode.setCompressOnPut(True)

        for propName in losLaserPropNames:
            for j in range(n_lasers):
                propNode_j = chanNode_i.addNode( \
                              laserName + '_%d_' % (j+1) + propName, 'numeric')
                propNode_j.setWriteOnce(False)
                propNode_j.setCompressOnPut(True)

        for propName in signalMappingPropNames:
            propNode = chanNode_i.addNode(propName, 'text')
            propNode.setWriteOnce(False)
            propNode.setCompressOnPut(True)

        channelNode = chanNode_i.addNode(realTimeChannelName, 'numeric')
        channelNode.setWriteOnce(False)
        channelNode.setCompressOnPut(True)

        polyNode = chanNode_i.addNode(polychromatorName, 'structure')
        for j in range(n_polychromators):
            polyNode_j = polyNode.addNode( \
                             polychromatorName + '_%d' % (j+1), 'structure')
            for propName in polychromatorPropNames:
                propNode = polyNode_j.addNode(propName, 'numeric')
                propNode.setWriteOnce(False)
                propNode.setCompressOnPut(True)
            for propName in polychromatorMappingPropNames:
                propNode = polyNode_j.addNode(propName, 'text')
                propNode.setWriteOnce(False)
                propNode.setCompressOnPut(True)

            # Line 568: adds nodes for empty set "POLYCHROMATORSIGNALSLASER"
            

if __name__ == '__main__':

    if len(sys.argv) != 2 or sys.argv[1] == '--help' or sys.argv[1] == '-h':
        print('Usage: python generate_nstxu_rtmpts_tree.py [tree_name]')
        print('Note: to create a tree with the name [tree_name], it is first ' \
               + 'necessary to ')
        print('    create an environment variable [tree_name]_path ' \
               + 'containing the location ')
        print('    where the tree files are to be stored.')
        exit()

    treeName = sys.argv[1]
    make_tree(treeName)


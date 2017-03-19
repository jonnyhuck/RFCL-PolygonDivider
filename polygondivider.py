# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PolygonDivider
								 A QGIS plugin
 Divides Polygons
							  -------------------
		begin				 : 2017-01-22
		git sha				 : $Format:%H$
		copyright			 : (C) 2017 by Roy Ferguson Consulting
		email				 : jonnyhuck@gmail.com
 ***************************************************************************/

/***************************************************************************
 *																		   *
 *	 This program is free software; you can redistribute it and/or modify  *
 *	 it under the terms of the GNU General Public License as published by  *
 *	 the Free Software Foundation; either version 2 of the License, or	   *
 *	 (at your option) any later version.								   *
 *																		   *
 ***************************************************************************/
"""
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication, Qt
from PyQt4.QtGui import QAction, QIcon, QFileDialog, QProgressBar
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from polygondivider_dialog import PolygonDividerDialog
import os.path
from squareishPolygonDividerQGIS import runSplit
from qgis.core import QgsVectorLayer, QgsMapLayerRegistry
from qgis.gui import QgsMessageBar

class PolygonDivider:
	"""QGIS Plugin Implementation."""
	
	progress = QProgressBar()

	def __init__(self, iface):
		"""Constructor.

		:param iface: An interface instance that will be passed to this class
			which provides the hook by which you can manipulate the QGIS
			application at run time.
		:type iface: QgsInterface
		"""
		# Save reference to the QGIS interface
		self.iface = iface
		# initialize plugin directory
		self.plugin_dir = os.path.dirname(__file__)
		# initialize locale
		locale = QSettings().value('locale/userLocale')[0:2]
		locale_path = os.path.join(
			self.plugin_dir,
			'i18n',
			'PolygonDivider_{}.qm'.format(locale))

		if os.path.exists(locale_path):
			self.translator = QTranslator()
			self.translator.load(locale_path)

			if qVersion() > '4.3.3':
				QCoreApplication.installTranslator(self.translator)


		# Declare instance attributes
		self.actions = []
		self.menu = self.tr(u'&Polygon Divider')
		# TODO: We are going to let the user set this up in a future iteration
		self.toolbar = self.iface.addToolBar(u'PolygonDivider')
		self.toolbar.setObjectName(u'PolygonDivider')
		

	# noinspection PyMethodMayBeStatic
	def tr(self, message):
		"""Get the translation for a string using Qt translation API.

		We implement this ourselves since we do not inherit QObject.

		:param message: String for translation.
		:type message: str, QString

		:returns: Translated version of message.
		:rtype: QString
		"""
		# noinspection PyTypeChecker,PyArgumentList,PyCallByClass
		return QCoreApplication.translate('PolygonDivider', message)


	def add_action(
		self,
		icon_path,
		text,
		callback,
		enabled_flag=True,
		add_to_menu=True,
		add_to_toolbar=True,
		status_tip=None,
		whats_this=None,
		parent=None):
		"""Add a toolbar icon to the toolbar.

		:param icon_path: Path to the icon for this action. Can be a resource
			path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
		:type icon_path: str

		:param text: Text that should be shown in menu items for this action.
		:type text: str

		:param callback: Function to be called when the action is triggered.
		:type callback: function

		:param enabled_flag: A flag indicating if the action should be enabled
			by default. Defaults to True.
		:type enabled_flag: bool

		:param add_to_menu: Flag indicating whether the action should also
			be added to the menu. Defaults to True.
		:type add_to_menu: bool

		:param add_to_toolbar: Flag indicating whether the action should also
			be added to the toolbar. Defaults to True.
		:type add_to_toolbar: bool

		:param status_tip: Optional text to show in a popup when mouse pointer
			hovers over the action.
		:type status_tip: str

		:param parent: Parent widget for the new action. Defaults None.
		:type parent: QWidget

		:param whats_this: Optional text to show in the status bar when the
			mouse pointer hovers over the action.

		:returns: The action that was created. Note that the action is also
			added to self.actions list.
		:rtype: QAction
		"""

		# Create the dialog (after translation) and keep reference
		self.dlg = PolygonDividerDialog()

		icon = QIcon(icon_path)
		action = QAction(icon, text, parent)
		action.triggered.connect(callback)
		action.setEnabled(enabled_flag)

		if status_tip is not None:
			action.setStatusTip(status_tip)

		if whats_this is not None:
			action.setWhatsThis(whats_this)

		if add_to_toolbar:
			self.toolbar.addAction(action)

		if add_to_menu:
			self.iface.addPluginToMenu(
				self.menu,
				action)

		self.actions.append(action)

		return action

	def initGui(self):
		"""Create the menu entries and toolbar icons inside the QGIS GUI."""

		icon_path = ':/plugins/PolygonDivider/icon.png'
		self.add_action(
			icon_path,
			text=self.tr(u'Divide Polygons'),
			callback=self.run,
			parent=self.iface.mainWindow())


	def unload(self):
		"""Removes the plugin menu item and icon from QGIS GUI."""
		for action in self.actions:
			self.iface.removePluginMenu(
				self.tr(u'&Polygon Divider'),
				action)
			self.iface.removeToolBarIcon(action)
		# remove the toolbar
		del self.toolbar


	def run(self):
		"""
		* Run method that performs all the real work
		"""

		# JJH: set up the dialog here--------------------------------------------
	
		# populate comboBox with the active layers
		self.dlg.comboBox.clear()	# need to clear here or it will add them all again every time the dialog is opened
		layers = self.iface.legendInterface().layers()
		layer_list = []
		for layer in layers:
			 layer_list.append(layer.name())
		self.dlg.comboBox.addItems(layer_list)

		# populate comboBox_2 with the possible directions
		self.dlg.comboBox_2.clear()	# need to clear here or it will add them all again every time the dialog is opened
		self.dlg.comboBox_2.addItems(['left to right', 'right to left', 'bottom to top', 'top to bottom'])

		# launch file browser for output file button - link to function
		self.dlg.pushButton.clicked.connect(self.select_output_file)

		#----------------------------------------------------------------------JJH

		# show the dialog
		self.dlg.show()
		
		# Run the dialog event loop
		result = self.dlg.exec_()
		
		# See if OK was pressed
		if result:

			# JH: RUN THE TOOL------------------------------------------------
			
			# add the progress bar			
			progressMessageBar = self.iface.messageBar().createMessage("Dividing Polygons...")
			self.progress.setMaximum(100)
			self.progress.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
			progressMessageBar.layout().addWidget(self.progress)
			self.iface.messageBar().pushWidget(progressMessageBar, level=QgsMessageBar.INFO)

			# get user settings
			inFile = layers[self.dlg.comboBox.currentIndex()]
			outFilePath = self.dlg.lineEdit_2.text()
			targetArea = float(self.dlg.lineEdit.text())
			absorbFlag = self.dlg.checkBox.isChecked()
			direction = self.dlg.comboBox_2.currentIndex()

			# run the tool
			runSplit(self, inFile, outFilePath, targetArea, absorbFlag, direction, self.progress)

			# add the result to the workspace
			layer = QgsVectorLayer(outFilePath, 'Divided Polygon', 'ogr')
			if layer.isValid():
				 QgsMapLayerRegistry.instance().addMapLayer(layer)	
			else:
				self.iface.messageBar().pushMessage("Error", "Failed to open resulting layer", level=QgsMessageBar.CRITICAL)
				
			# reset progress bar
			self.dlg.progressBar.setValue(0)
			
			# remove progress bar from the message bar and add in a success message
			self.iface.messageBar().clearWidgets()
			self.iface.messageBar().pushMessage("Success!", "Polygon Divided Successfully!", level=QgsMessageBar.INFO)

			#--------------------------------------------------------------JJH


	def select_output_file(self):
		"""
		* JJH:Open file browser
		"""

		# get filename from dialog	
		filename = QFileDialog.getSaveFileName(self.dlg, "Select output file ","", '*.shp')
	
		# verify that a name was selected
		if filename != "":

			# clear previous value
			self.dlg.lineEdit_2.clear()

			# make sure that an extension was included
			if filename[-4:] != '.shp':
				filename += '.shp'

			# put the result in the text box on the dialog
			self.dlg.lineEdit_2.setText(filename)
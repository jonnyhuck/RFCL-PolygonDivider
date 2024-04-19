from qgis import processing
from processing.gui.AlgorithmExecutor import execute_in_place
from qgis.core import QgsApplication


def rotate_layer_in_place(input_layer, degrees: float, rotation_center: str = None):
    """rotates layer in place"""
    registry = QgsApplication.instance().processingRegistry()
    alg = registry.algorithmById("native:rotatefeatures")

    params = {
        'INPUT': input_layer,
        'ANGLE': degrees,
        'ANCHOR': rotation_center
    }

    execute_in_place(alg, params)


def rotate_layer(input_path: str, degrees: float, rotation_center: str = None, output_path: str = 'TEMPORARY_OUTPUT'):
    """returns rotated in-memory layer"""
    return processing.run("native:rotatefeatures",
        {
            'INPUT': input_path,
            'ANGLE': degrees,
            'ANCHOR': rotation_center,
            'OUTPUT': output_path
        }
    )['OUTPUT']

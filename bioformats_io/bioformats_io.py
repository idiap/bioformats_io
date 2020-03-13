"""
- Bioformat loader and converter from TIFF, BigTIFF, OME-TIFF, STK, LSM, LIF, SGI, ImageJ, MicroManager,
FluoView, SEQ, GEL, ... files to OME-TIFF specifications with metadatas.

- Context manager for the calling/closing of a jvm when using
bioformats

Copyright (c) 2008 Idiap Research Institute, http://www.idiap.ch/
Written by Olivia Mariani <olivia.mariani@idiap.ch>
Adrian Shajkofci <adrian.shajkofci@idiap.ch>

This file is part of Unshearing.

Unshearing is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

Unshearing is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Foobar. If not, see <http://www.gnu.org/licenses/>.

"""

import numpy as np
import os
import logging
import argparse

import javabridge
import bioformats
import bioformats.omexml as ome


def open_loci_5d(input_file_path, output_file_path='', nt=-1, nc=-1, nz=-1, ti=0, ci=0, zi=0, series=0,
                 logging_level='INFO', bioformats_logging='OFF'):
    """
    Bioformats file reader.
    Warning: The data is expected to be in XYZCT order.
    :param input_file_path: Input file path (OME, OMETIFF, CZI, LIF, TIFF, LIF, ...)
    :param output_file_path: Output NPY file path. By default: input_file_path with NPY extension.
    :param nt: Number of time points
    :param nc: Number of channels
    :param nz: Number of slices
    :param ti: Initial time point
    :param ci: Initial channel
    :param zi: Initial slice
    :param logging_level: Level of logging: DEBUG, INFO, WARNING, ERROR
    :param bioformats_logging: Bioformats logging level
    :param series: Which series to open, if more than one
    :return: 5D ndarray in order XYZCT
    """
    dir_data = os.path.dirname(input_file_path)
    logging.basicConfig(filename=os.path.join(dir_data, 'loci_reader.log'), level=logging_level)

    logging.info('Opening file {}...'.format(input_file_path))
    loci_logging(bioformats_logging)
    reader = bioformats.get_image_reader(input_file_path, path=input_file_path)
    image_reader = reader.rdr

    nx = image_reader.getSizeX()
    ny = image_reader.getSizeY()

    nz_tot = image_reader.getSizeZ()
    nc_tot = image_reader.getSizeC()
    nt_tot = image_reader.getSizeT()

    if nz == -1:
        nz = nz_tot-zi
    if nc == -1:
        nc = nc_tot - ci
    if nt == -1:
        nt = nt_tot - ti

    if not (0 <= zi < nz_tot):
        logging.warning('The initial slice {} is not within the slices range 0-{}. Using central slice {} instead.'
                        .format(zi, nz_tot-1, int(nz_tot*0.5)))
        zi = int(nz_tot*0.5)

    if nz_tot < nz+zi:
        logging.warning('Total number of slices is {} and you chose {}. Using {} instead.'
                        .format(nz_tot, nz+zi, nz_tot-zi))
        nz = nz_tot-zi
        
    if not (0 <= ci < nc_tot):
        logging.warning('The initial channel {} is not within the channel range 0-{}. Using central slice {} instead.'
                        .format(ci, nc_tot - 1, int(nc_tot * 0.5)))
        ci = int(nc_tot * 0.5)
    if nc_tot < nc + ci:
        logging.warning('Total number of channel is {} and you chose {}. Using {} instead.'
                        .format(nc_tot, nc + ci, nc_tot - ci))
        nc = nc_tot - ci
               
    if not (0 <= ti < nt_tot):
        logging.warning('The initial time point {} is not within the time-points range {}-{}. '
                        'Using time point 0 instead.'.format(ti, 0, nt_tot-1))
        ti = 0
    if nt_tot < nt+ti:
        logging.warning('Total number of time points is {} and you chose {}. Using {} instead.'
                        .format(nt_tot, nt+ti, nt_tot-ti))
        nt = nt_tot-ti

    if not (0 <= series < image_reader.getSeriesCount()):
        logging.warning('The series {} is not within the series scope 0-{}. Opening series 0 instead.'
                        .format(series, image_reader.getSeriesCount()))
        series = 0

    image_reader.getSeries(series)
    format_tool = bioformats.formatreader.make_format_tools_class()
    pixel_type = format_tool.getPixelTypeString(reader.rdr.getPixelType())
    image5d = np.zeros([ny, nx, nz - zi, nc - ci, nt - ti], pixel_type)

    for t in np.arange(ti, nt+ti):
        logging.debug('Opening time point {} / {}'.format(t, nt))
        for c in np.arange(ci, nc+ci):
            for z in np.arange(zi, nz+zi):
                image5d[..., z - zi, c - ci, t - ti] = reader.read(c=c, z=z, t=t, series=series, rescale=False)

    dir_data = os.path.dirname(input_file_path)
    filename = os.path.basename(input_file_path)

    if output_file_path == '':
        output_file_path = os.path.join(dir_data, filename[:filename.find(".")] + '.npy')

    np.save(output_file_path, image5d)
    logging.info('Successfully save file {}'.format(output_file_path))


def save_to_ome(input_file_path, output_file_path='', file_name_for_metadata=None, resolutions=np.array([]),
                units=np.array([]), pixel_type='float', logging_level='INFO', bioformats_logging='OFF'):
    """
    :param input_file_path: Input NPY file. Expected XYZCT 5-D data.
    :param output_file_path: /path/to/output.ome Can be OME or OMETIFF. By default: Input file name with changed extension.
    :param file_name_for_metadata: Will copy the metadata info of this file to the output file
    :param resolutions: Pixel physical size can be inserted here
    :param units: Pixel physical size units
    :param pixel_type: Pixel type
    :param logging_level: level of logging: DEBUG, INFO, WARNING, ERROR
    :param bioformats_logging: Bioformats logging level
    :return:
    """
    dir_data = os.path.dirname(input_file_path)
    filename = os.path.basename(input_file_path)

    logging.basicConfig(filename=os.path.join(dir_data, 'loci_writer.log'), level=logging_level)

    if output_file_path == '':
        output_file_path = os.path.join(dir_data, filename[:filename.find(".")]+'.ome')

    logging.info('Saving to {} ...'.format(output_file_path))

    loci_logging(bioformats_logging)

    image5d = np.load(input_file_path)

    image_writer = bioformats.formatwriter.make_image_writer_class()
    writer = image_writer()

    if file_name_for_metadata is not None:
        rdr = bioformats.ImageReader(file_name_for_metadata, perform_init=True)
        jmd = javabridge.JWrapper(rdr.rdr.getMetadataStore())
        pixel_type = jmd.getPixelsType(0).getValue()

        omexml = ome.OMEXML()
        omexml.image(0).Name = os.path.split(output_file_path)[1]
        p = omexml.image(0).Pixels
        assert isinstance(p, ome.OMEXML.Pixels)

        p.node.set("PhysicalSizeX", str(jmd.getPixelsPhysicalSizeX(0).value()))
        p.node.set("PhysicalSizeY", str(jmd.getPixelsPhysicalSizeY(0).value()))
        p.node.set("PhysicalSizeZ", str(jmd.getPixelsPhysicalSizeZ(0).value()))
        p.node.set("PhysicalSizeXUnit", jmd.getPixelsPhysicalSizeX(0).unit().getSymbol())
        p.node.set("PhysicalSizeYUnit", jmd.getPixelsPhysicalSizeY(0).unit().getSymbol())
        p.node.set("PhysicalSizeZUnit", jmd.getPixelsPhysicalSizeZ(0).unit().getSymbol())

        p.SizeX = image5d.shape[1]
        p.SizeY = image5d.shape[0]
        p.SizeC = image5d.shape[3]
        p.SizeT = image5d.shape[4]
        p.SizeZ = image5d.shape[2]

        p.PixelType = pixel_type
    else:
        resolutions = np.array(resolutions)
        units = np.array(units)
        omexml = ome.OMEXML()
        omexml.image(0).Name = os.path.split(output_file_path)[1]
        p = omexml.image(0).Pixels
        assert isinstance(p, ome.OMEXML.Pixels)
        if resolutions.size:
            rx, ry, rz = resolutions
            p.node.set("PhysicalSizeX", str(rx))
            p.node.set("PhysicalSizeY", str(ry))
            p.node.set("PhysicalSizeZ", str(rz))
        if units.size:
            ux, uy, uz = units
            p.node.set("PhysicalSizeXUnit", ux)
            p.node.set("PhysicalSizeYUnit", uy)
            p.node.set("PhysicalSizeZUnit", uz)
        p.SizeX = image5d.shape[1]
        p.SizeY = image5d.shape[0]
        p.SizeC = image5d.shape[3]
        p.SizeZ = image5d.shape[2]
        p.SizeT = image5d.shape[4]

        p.PixelType = pixel_type
        p.channel_count = image5d.shape[3]

    xml = omexml.to_xml()
    script = """
    importClass(Packages.loci.formats.services.OMEXMLService,
                Packages.loci.common.services.ServiceFactory,
                Packages.loci.formats.MetadataTools,
                Packages.loci.formats.meta.IMetadata);

    service = new ServiceFactory().getInstance(OMEXMLService);
    metadata = service.createOMEXMLMetadata(xml);

    """
    meta = javabridge.run_script(script, dict(xml=xml))
    if os.path.isfile(output_file_path):
        os.remove(output_file_path)

    writer.setMetadataRetrieve(meta)
    writer.setId(output_file_path)

    nz, nc, nt = image5d.shape[2:]

    for t in range(nt):
        logging.debug('Saving frame {} ...'.format(t))
        for c in range(nc):
            for z in range(nz):
                index = z + nz * c + nz * nc * t
                save_im = bioformats.formatwriter.convert_pixels_to_buffer(image5d[..., z, c, t], p.PixelType)
                writer.saveBytesIB(index, save_im)
    writer.close()
    logging.info('File {} successfully saved.'.format(output_file_path))


def loci_logging(logging_level='OFF'):
    """
    Call bioformats library's logger Apache Log4j and disable debug logging (and more specifically debug checks,
    which are very slow)
    :param logging_level: can be 'ALL', 'DEBUG', 'ERROR', 'FATAL', 'INFO', 'OFF', 'TRACE', 'WARN'
    """

    script = """
    importClass(Packages.loci.common.DebugTools);
    DebugTools.enableLogging(logging_level);
    """

    javabridge.run_script(script, bindings_in={"logging_level": logging_level})


def parsing():
    """
    Bash commands. Mendatory: Input and Output absolute file path
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file_path", type=str,
                        help="Input data path. Expected 5D data XYZCT.")
    parser.add_argument("-o", "--output_file_path", type=str, default='',
                        help="Output data path.")
    parser.add_argument("-z", "--z", type=int, default="0",
                        help="Slice, if more than one in the input file.")
    parser.add_argument("-c", "--c", type=int, default="0",
                        help="Channel, if more than one in the input file.")
    parser.add_argument("--ti", type=int, default="0",
                        help="Initial time point. Default=0")
    parser.add_argument("--nt", type=int, default="-1",
                        help="Number of time points. Default=-1, opening from -ti until last time point.")
    parser.add_argument("--series", type=int, default="0",
                        help="Series, if more than one in the input file.")
    parser.add_argument("--heap_size", type=str, default="8G",
                        help="Maximum Java heap size. Default=8G.")
    parser.add_argument("--save_ome", type=bool, default=False,
                        help="Saves the data as OME files. Saved by default in input_file_path directory.")
    parser.add_argument("--metadata_path", type=str,
                        help="Path to the original bioformats file, to retrieves the metadata.")
    parser.add_argument("--xyz_resolutions_float", type=str, default="0.1,0.1,0.1", nargs=3,
                        help="x,y, and z resolutions. Example: --xyz_resolutions_float 1.56 1.56 3.56")
    parser.add_argument("--xyz_unit_str", type=str, default="µm,µm,µm",
                        help="x,y, and z units. Example: --xyz_resolutions_float µm,µm,µm")
    parser.add_argument("--pixel_type", type=str, default='float',
                        help="Pixel type. Default=float")
    parser.add_argument("--logging", type=str,
                        default='INFO',
                        help="Logging levels. Can be: INFO, WARNING, or 'ERROR'. Default: INFO")
    parser.add_argument("--bioformats_logging", type=str, default='OFF',
                        help="Bioformats logging. Can be ALL', 'DEBUG', 'ERROR', 'FATAL', 'INFO', 'OFF', 'TRACE', "
                             "'WARN'")
    return parser.parse_args()


class JVM:

    def __init__(self, java_heap_size):
        self.java_heap_size = java_heap_size

    def __enter__(self):
        javabridge.start_vm(class_path=bioformats.JARS, max_heap_size=self.java_heap_size)

    def __exit__(self, *args):
        javabridge.kill_vm()


if __name__ == "__main__":

    parse = parsing()

    with JVM(parse.heap_size):

        logging.info('')

        if parse.save_ome:

            res = np.array([float(item)for item in parse.xyz_resolutions_float.split(',')])
            un = np.array([str(item) for item in parse.xyz_unit_str.split(',')])

            save_to_ome(parse.input_file_path, output_file_path=parse.output_file_path,
                        file_name_for_metadata=parse.metadata_path, resolutions=res, units=un,
                        pixel_type=parse.pixel_type, logging_level=parse.logging,
                        bioformats_logging=parse.bioformats_logging)
        else:
            open_loci_5d(parse.input_file_path, output_file_path=parse.output_file_path, nt=parse.nt, nc=1, nz=1,
                         ti=parse.ti, ci=parse.c, zi=parse.z, series=parse.series, logging_level=parse.logging,
                         bioformats_logging=parse.bioformats_logging)

#!/usr/bin/env python
#

"""
Classes and utilities that define a set of regions of interest to tile the sky
"""

import sys
import os

import healpy
import numpy as np
import yaml

from astropy import table
from astropy import coordinates

from fermipy import skymap
from fermipy import wcs_utils


def extract_map_data_from_file(fitsfile, lons, lats, **kwargs):
    """ Extract values for a set of coordinates from a fit file

    Parameters
    ----------
    lons : array-like
        'longitudes' (RA or GLON) in degrees

    lats : array-like
        'latitudes' (DEC or GLAT) in degrees

    Returns
    ----------
    vals : numpy.ndarray((n)) or numpy.ndarray((n,nbins))
           Values of pixels in the flattened map, np.nan used to flag coords outside of map
           Depending on the shape of the FITS HDU vals could be either 1D or 2D.
    """
    try:
        mdata = skymap.Map.create_from_fits(fitsfile, **kwargs)
    except OSError:
        print("Failed to read %s" % fitsfile)
        nbin = kwargs.get('nbin', None)
        if nbin is None:
            return np.zeros((len(lons)))
        else:
            return np.zeros((len(lons), nbin))

    ibin = kwargs.get('ibin', None)
    vals = mdata.get_map_values(lons, lats, ibin)

    return vals


def extract_table_data_from_file(fitsfile, lons, lats, **kwargs):
    """ Extract values for a set of coordinates from a fit file

    Parameters
    ----------
    lons : array-like
        'longitudes' (RA or GLON) in degrees

    lats : array-like
        'latitudes' (DEC or GLAT) in degrees

    Returns
    -------
    vals : numpy.ndarray((n)) or numpy.ndarray((n,nbins))
           Values of pixels in the flattened map, np.nan used to flag coords outside of map
           Depending on the shape of the FITS HDU vals could be either 1D or 2D.
    """
    tablelist = kwargs.get('tables', None)
    try:
        mdata = skymap.Map.create_from_fits(fitsfile, **kwargs)
        if tablelist is None:
            tab = table.Table.read(fitsfile)
        else:
            tlist = []
            for tname in tablelist:
                tlist.append(table.Table.read(fitsfile, tname))
            tab = table.hstack(tlist)
    except OSError:
        print("Failed to read %s" % fitsfile)
        return None
    pix_idxs = mdata.get_pixel_indices(lons, lats)
    out_table = tab[pix_idxs]
    return out_table


def extract_column_data_from_file(fitsfile, lons, lats, **kwargs):
    """ Extract values for a set of coordinates from a fit file

    Parameters
    ----------
    lons : array-like
        'longitudes' (RA or GLON) in degrees

    lats : array-like
        'latitudes' (DEC or GLAT) in degrees

    tabname : str
        name of the table containing the column

    colname : str
        name of the column

    Returns
    -------
    vals : numpy.ndarray((n)) or numpy.ndarray((n,nbins))
           Values of pixels in the flattened map, np.nan used to flag coords outside of map
           Depending on the shape of the FITS HDU vals could be either 1D or 2D.
    """
    tabname = kwargs.get('tabname', None)
    colname = kwargs.get('colname', None)
    try:
        mdata = skymap.Map.create_from_fits(fitsfile, **kwargs)
        if tabname is None:
            tab = table.Table.read(fitsfile)
        else:
            tab = table.Table.read(fitsfile, tabname)
    except OSError:
        nbin = kwargs.get('nbin', None)
        if nbin is None:
            return np.zeros((len(lons)))
        else:
            return np.zeros((len(lons), nbin))
    pix_idxs = mdata.get_pixel_indices(lons, lats)
    out_table = tab[colname][pix_idxs]
    return out_table


def merge_and_sort_table(tdict, klist):
    """ Build a single table by pulling a list of rows out of a dict of tables
    """
    ilocal = 0
    kdict = {}
    rlist = []

    # Sort the rows
    for key in klist:
        if kdict.has_key(key):
            kdict[key] += 1
        else:
            kdict[key] = 0
        ilocal = kdict[key]
        rlist.append(tdict[key][ilocal])

    if len(rlist) == 0:
        return None

    # This builds the columns from the first entry, and removes all the rows
    out_table = table.Table(rlist[0].columns)
    to_remove = range(len(out_table))
    out_table.remove_rows(to_remove)

    # Now loop on the rList and add back the rows
    for row in rlist:
        out_table.add_row(row)
    return out_table


def circum_at_lat(lat):
    """ Get the circumference of a ring at a given latitude

    Parameters
    ----------
    lat : array-like
        'latitude' (DEC or GLAT) in radians

    Returns
    -------
    Circumference (in radians)
    """
    return 2. * np.pi * np.cos(lat)


def max_circums(lat_list):
    """ Get the maximum circumference for each set of latitudes

    Parameters
    ----------
    lat_list : list of array-like objects
       sets of 'latitude' (DEC or GLAT) in radians

    Returns
    -------
    circ_array : array-like
       Maximum circumference (in radians)from each set
    """
    lat_array = np.vstack([lat_list])
    circ_array = circum_at_lat(lat_array)
    return circ_array.max(0)


def ring_lats(nrings):
    """ Get a set of evenly spaced latitudes

    Parameters
    ----------
    nrings :  int
       number of rings

    Returns
    -------
    lat_edges : np.ndarray((nrings+1))
       edges of rings (in radians)
    """
    lat_edges = np.linspace(-0.5 * np.pi, 0.5 * np.pi, nrings + 1)
    return lat_edges


def ring_lats_from_roi(roisize):
    """  Get a set of evenly spaced latitudes

    Parameters
    ----------
    roisize : float
       minimum size of region of interest (in radians)

    Returns
    -------
    lat_edges : np.ndarray((nrings+1))
       edges of rings (in radians)
    """
    nrings = int(np.ceil(np.pi / roisize))
    return ring_lats(nrings)


def ring_lons(max_circ, roisize):
    """ Get a set of evenly spaced longitudes

    Parameters
    ----------
    max_circs : float
       Maximum circumeference of the ring

    roisize : float
       minimum size of region of interest (in radians)

    Returns
    -------
    lon_edges : np.ndarray((nroi+1))
       edges of rois (in radians)
    """
    nrings = int(np.ceil(max_circ / roisize))
    lon_edges = np.linspace(-1.0 * np.pi, np.pi, nrings + 1)
    return lon_edges


def get_roi_centers(roisize):
    """ Get the centers of a set of rois that are evenly spaced

    Parameters
    ----------
    roisize : float
       minimum size of region of interest (in radians)

    Returns
    -------
    roi_lats : np.array (nroi)
       'latitudes' of centers of ROIs

    roi_lons : np.array (nroi)
       'longitudes' of centers of ROIs

    n_roi_per_ring : np.ndarray((n_rings),'i')
       number of ROI in each ring
    """
    lat_edges = ring_lats_from_roi(roisize)
    lat_cents = (lat_edges[0:-1] + lat_edges[1:]) / 2.
    n_rings = len(lat_cents)
    max_circs = max_circums([lat_edges[0:-1], lat_cents, lat_edges[1:]])

    lat_cent_list = []
    lon_cent_list = []

    n_roi_per_ring = np.ndarray((n_rings), 'i')
    for idx, (lat_cent, max_circ) in enumerate(zip(lat_cents, max_circs)):
        lon_edges = ring_lons(max_circ, roisize)
        lon_cents = (lon_edges[0:-1] + lon_edges[1:]) / 2.
        n_roi_per_ring[idx] = len(lon_cents)
        for lon_cent in lon_cents:
            lat_cent_list.append(lat_cent)
            lon_cent_list.append(lon_cent)

    roi_lats = np.array(lat_cent_list)
    roi_lons = np.array(lon_cent_list)

    return roi_lats, roi_lons, n_roi_per_ring


def lat_to_rings(roisize_deg, lats):
    """ Get the rings corresponding to a set of latitidues

    Parameters
    ----------
    roisize_def : float
       minimum size of region of interest (in degress)

    lats : array-like
       input set of latitudes ( in degress )

    Returns
    -------
    rings : np.ndarray((n_rings),'i')
       Index of corresponding ring for each latitude
    """
    nrings = int(np.ceil(180. / roisize_deg))
    roi_size = 180. / nrings
    rings = (np.floor((lats + 90.) / roi_size)).astype(int)
    return rings


def lon_idxs(rings, lons, nroi_per_ring):
    """ Get the ROI indices corresponding to a set of rings and longitudes

    Parameters
    ----------
    rings : array-like
       Indices of corresponding rings

    lons : array-like
       input set of longitudes ( in degrees )

    Returns
    -------
    idxs : np.ndarray((n_roi),'i')
       Indices of corresponding rois
    """
    nrois = nroi_per_ring[rings]
    ring_idxs = np.zeros((len(nroi_per_ring)), 'i')
    ring_idxs[1:] = nroi_per_ring.cumsum()[0:-1]
    roisizes = 360. / nrois

    lons_fold = np.where(lons > 180., lons - 180., lons + 180)
    roi_idxs = (np.floor((lons_fold) / roisizes)).astype(int) + ring_idxs[rings]
    return roi_idxs


class DMROISet(object):
    """ A simple class to generate the centers of a set of ROIs tile the entire sky

    This is just the base class.  Implementations of the tilings are in sub-classes:

    DMROISetHPX : HEALPix-based tiling
    DMROISetWCS : Cartesian tiling
    """

    def __init__(self, prefix='roi_'):
        """ C'tor, takes the prefix used to construct the ROI names
        """
        self.__prefix = prefix
        self.__roi_dict = None
        self._nroi = 0
        self._glats = None
        self._glons = None

    @property
    def glats(self):
        """ The set of Galactic latitudes
        """
        return self._glats

    @property
    def glons(self):
        """ The set of Galactic longitudes
        """
        return self._glons

    @property
    def nroi(self):
        """ The number of ROIs
        """
        return self._nroi

    @property
    def roi_dict(self):
        """ A dictionary of dictionaries, each of which can be used to create an ROI
        """
        return self.__roi_dict

    @property
    def prefix(self):
        """ The prefix used in ROI name construction
        """
        return self.__prefix

    @staticmethod
    def create_from_yaml(yamlfile):
        """ Create a DMROISet from a yaml file
        """
        din = yaml.safe_load(open(yamlfile))
        prefix = din['prefix']
        basedir = din['basedir']
        if din['proj_type'].lower() in ['healpix', 'hpx']:
            nside = din['hpx_nside']
            nest = din['hpx_nest']
            roiset_out = DMROISetHPX(nside, nest, prefix)
        elif din['proj_type'].lower() in ['wcs']:
            roisize_deg = din['roi_size']
            roiset_out = DMROISetWCS(roisize_deg, prefix)
        return roiset_out, basedir

    def _build_roi_dict(self):
        """ Build the dictionary of dictionaries from set of latitudes and longitudes
            provided by the sub-class.
        """
        self.__roi_dict = {}
        self._nroi = len(self._glats)
        roi_idxs = np.arange(self._nroi)
        for glat, glon, roi_idx in zip(self._glats, self._glons, roi_idxs):
            self.__roi_dict["%s%05i" % (self.__prefix, roi_idx)] = (
                {"selection": {"glat": float(glat), "glon": float(glon)}})
        return self.__roi_dict

    def write_yaml(self, filepath):
        """ Write the dictionary of dictionaries into a yaml file
        """
        fout = open(filepath, "w")
        yaml.dump(self.__roi_dict, fout)
        fout.close()

    def lb_to_idx(self, glons, glats):
        """ Convert from Galactic coordinates to ROI index
        """
        raise NotImplementedError("DMROISet.lb_to_idx")


class DMROISetWCS(DMROISet):
    """ A sub-class of DMROISet that uses a simple cartensian grid tiling
    """

    def __init__(self, roisize_deg, prefix='roi_'):
        """ C'tor, specify the size of each ROI in degrees
        """
        DMROISet.__init__(self, prefix)
        self.__roisize_deg = roisize_deg

        glats_rad, glons_rad, self.__nroi_per_ring = get_roi_centers(np.radians(self.__roisize_deg))
        self.__ring_idx = np.zeros((len(self.__nroi_per_ring)), 'i')
        self.__ring_idx[1:] = self.__nroi_per_ring.cumsum()[0:-1]
        self._glats = np.degrees(glats_rad)
        self._glons = np.degrees(glons_rad)

        self._build_roi_dict()

    @property
    def nroi_per_ring(self):
        """ Return the number of ROIs for each ring """
        return self.__nroi_per_ring

    @property
    def ring_idx(self):
        """ Return the indices for the first ROI in each ring """
        return self.__ring_idx

    @property
    def roisize_deg(self):
        """ The size of each roi, in degrees

        The ROIs are square, and this gives the length of a side
        """
        return self.__roisize_deg

    def lb_to_idx(self, glons, glats):
        """ Get the ROI indices corresponding to a set of rings and longitudes

        Parameters
        ----------
        glons : array-like
           Galactic longitudes

        glats : array-like
           Galactic latitdues

        Returns
        -------
        idxs : np.ndarray((n_roi),'i')
            Indices of corresponding rois
        """
        rings = lat_to_rings(self.__roisize_deg, glats)
        roi_idxs = lon_idxs(rings, glons, self.__nroi_per_ring)
        return roi_idxs

    def create_roi_wcs(self, roi_idx, coordsys='CEL', projection='AIT',
                       cdelt=1.0, crpix=1., naxis=2, energies=None):
        """ Create a WCS object.

        Parameters
        ----------
        roi_idx : int
            Index of the corresponding ROI

        coordsys : str

        projection : str

        cdelt : float

        crpix : float or (float,float)
            In the first case the same value is used for x and y axes

        Returns
        -------
        wcs : `~astropy.wcs`
        """
        skydir = coordinates.SkyCoord(self._glons[roi_idx], self._glats[roi_idx],
                                      frame=coordinates.Galactic, unit="deg")
        wcs = wcs_utils.create_wcs(skydir, coordsys, projection, cdelt, crpix, naxis, energies)
        return wcs

    def extract_map_data(self, shape, sky_crds, filestr, basedir=".", **kwargs):
        """ Extract data from files built as an roi_set

        The data are assumed to be in a set of files with names like
           roi_00000/tsmap.fits

        Parameters
        ----------
        shape : tuple(int)
            shape of the output array.
            Either shape[0] or shape[0] x shape[1] must be equal to the number of directions

        sky_crds : array-like (2,ndir)
            Directions to extract data for

        filestr : str
            Name of file within each ROI sub-directory

        basedir : str
            Path to top-level directory

        Returns
        -------
        out_array : ~numpy.ndarray(shape)
            Array with the extracted data
        """
        nbin = kwargs.get('nbin', None)

        # Mask out the stuff outide the projections
        mask = np.invert(np.isnan(sky_crds)).any(1)
        good_sky = sky_crds[mask]
        if nbin is not None:
            full_mask = np.expand_dims(mask, -1) * np.ones((nbin), bool)
            full_mask = full_mask.reshape((full_mask.size))
        else:
            full_mask = mask

        # This is the array we are filling
        out_array = np.zeros(shape)
        out_array.flat[np.invert(full_mask.flat)] = np.nan

        # Figure out which ROI each pixel belongs to
        rois = self.lb_to_idx(good_sky[0:, 0], good_sky[0:, 1])

        # This is an array of only the pixels inside the projection
        if nbin is None:
            copy_array = np.zeros((mask.sum()))
        else:
            copy_array = np.zeros((mask.sum(), nbin))

        # Loop over ROIs
        min_roi = rois.min()
        max_roi = rois.max()
        for iroi in xrange(min_roi, max_roi + 1):
            # progress
            if iroi % 10 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()

            # Mask out the directions for this ROI
            local_mask = rois == iroi

            if not local_mask.any():
                continue

            # Open the corresonding file
            filepath = os.path.join(basedir, "%s%05i" % (self.prefix, iroi), filestr)

            # Extract the data for this ROI
            copy_array[local_mask] = extract_map_data_from_file(filepath,
                                                                good_sky[0:, 0][local_mask],
                                                                good_sky[0:, 1][local_mask],
                                                                **kwargs)

        # Copy the data into the full map
        out_array.flat[full_mask] = copy_array.flat
        return out_array

    def extract_table_data(self, sky_crds, filestr, basedir=".", **kwargs):
        """ Extract data from files built as an roi_set

        The data are assumed to be in a set of files with names like
           roi_00000/tsmap.fits

        Parameters
        ----------
        sky_crds : array-like (2,ndir)
            Directions to extract data for

        filestr : str
            Name of file within each ROI sub-directory

        basedir : str
            Path to top-level directory

        tables : [str...]
            List of names of the tables to extract

        Returns
        -------
        out_table : ~astropy.table.Table
            Table with the extracted data
        """
        # This is the table
        out_table = None

        # Figure out which ROI each direction belongs to
        rois = self.lb_to_idx(sky_crds[0], sky_crds[1])

        table_dict = {}

        # Loop over ROIs
        min_roi = rois.min()
        max_roi = rois.max()

        nrows = 0
        for iroi in xrange(min_roi, max_roi + 1):
            # progress
            if iroi % 10 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()

            # Mask out the directions for this ROI
            local_mask = rois == iroi

            if not local_mask.any():
                continue

            # Open the corresonding file
            filepath = os.path.join(basedir, "%s%05i" % (self.prefix, iroi), filestr)

            # Extract the data for this ROI
            table_data = extract_table_data_from_file(filepath,
                                                      sky_crds[0][local_mask],
                                                      sky_crds[1][local_mask],
                                                      **kwargs)
            table_dict[iroi] = table_data
            nrows += len(table_data)

        out_table = merge_and_sort_table(table_dict, rois)
        return out_table

    def extract_column_from_tables(self, shape, sky_crds, filestr, basedir=".", **kwargs):
        """
        Extract data from files built as an roi_set

        The data are assumed to be in a set of files with names like
           roi_00000/tsmap.fits

        Parameters
        ----------
        sky_crds : array-like (2,ndir)
            Directions to extract data for

        filestr : str
            Name of file within each ROI sub-directory

        basedir : str
            Path to top-level directory

        tabname : str
            Name of table containing data to extract

        colname : std
            Name of column to extract

        Returns
        -------
        out_array : ~numpy.ndarray(shape)
            Array with the extracted data
        """
        nbin = kwargs.get('nbin', None)

        # Mask out the stuff outide the projections
        mask = np.invert(np.isnan(sky_crds)).any(1)
        good_sky = sky_crds[mask]

        if nbin is not None:
            full_mask = np.expand_dims(mask, -1) * np.ones((nbin), bool)
            full_mask = full_mask.reshape((full_mask.size))
        else:
            full_mask = mask

        # This is the array we are filling
        out_array = np.zeros(shape)
        out_array.flat[np.invert(full_mask.flat)] = np.nan

        # Figure out which ROI each direction belongs to
        rois = self.lb_to_idx(good_sky[0:, 0], good_sky[0:, 1])

        # This is an array of only the pixels inside the projection
        if nbin is None:
            copy_array = np.zeros((mask.sum()))
        else:
            copy_array = np.zeros((mask.sum(), nbin))

        # Loop over ROIs
        min_roi = rois.min()
        max_roi = rois.max()

        for iroi in xrange(min_roi, max_roi + 1):
            # progress
            if iroi % 10 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()

            # Mask out the directions for this ROI
            local_mask = rois == iroi

            if not local_mask.any():
                continue

            # Open the corresonding file
            filepath = os.path.join(basedir, "%s%05i" % (self.prefix, iroi), filestr)

            # Extract the data for this ROI
            copy_array[local_mask] = extract_column_data_from_file(filepath,
                                                                   good_sky[0:, 0][local_mask],
                                                                   good_sky[0:, 1][local_mask],
                                                                   **kwargs)

        out_array.flat[full_mask] = copy_array.flat
        return out_array

    def extract_single_table(self, filestr, basedir=".", **kwargs):
        """ Extract data from a single table in the ROI set

        The data are assumed to be in a set of files with names like
           roi_00000/tsmap.fits

        Parameters
        ----------
        filestr : str
            Name of file within each ROI sub-directory

        basedir : str
            Path to top-level directory

        Returns
        -------
        out_table : ~astropy.table.Table
            Table with the extracted data
        """
        # Open the corresonding file
        iroi = 0
        filepath = os.path.join(basedir, "%s%05i" % (self.prefix, iroi), filestr)
        hdu = kwargs.get('table', 'EBOUNDS')
        tab = table.Table.read(filepath, hdu)
        return tab


class DMROISetHPX(DMROISet):
    """ A sub-class of DMROISet that uses a HEALPix tiling
    """

    def __init__(self, nside, nest=False, prefix='roi_'):
        """ C'tor, specify the HEALPix parameters
        """
        DMROISet.__init__(self, prefix)
        self.__nside = nside
        self.__nest = nest
        self.__nroi = 12 * nside * nside

        thetas, phis = healpy.pixelfunc.pix2ang(self.__nside, np.arange(self.__nroi), self.__nest)
        self._glats = np.degrees(0.5 * np.pi - thetas)
        self._glons = np.degrees(phis)

        self._build_roi_dict()

    @property
    def nside(self):
        """ The HEALPix nside parameter
        """
        return self.__nside

    @property
    def nest(self):
        """ The HEALPix parameter for NESTED or RING indexing
        """
        return self.__nest

    def lb_to_idx(self, glons, glats):
        """ Convert from Galactic coordinates to ROI index
        """
        raise NotImplementedError("DMROISetHPX.lb_to_idx")


if __name__ == '__main__':

    ROI_SIZE_DEG = 8.0  # in deg
    ROI_SET = DMROISetWCS(ROI_SIZE_DEG)
    ROI_IDXS = ROI_SET.lb_to_idx(ROI_SET.glons, ROI_SET.glats)

    # ROI_SET.write_yaml("rois.yaml")

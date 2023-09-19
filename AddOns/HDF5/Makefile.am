pkglib_LTLIBRARIES    =	\
	libSherpaHDF5v1Output.la \
	libSherpaHDF5Output.la

GITTAG = AddOns/HDF5
include ../../git.make

libSherpaHDF5v1Output_la_SOURCES = \
	Git_Info.C \
	Output_HDF5_v1.C

libSherpaHDF5Output_la_SOURCES = \
	Git_Info.C \
	Output_HDF5.C

localinc_HEADERS = \
	highfive/H5DataSet.hpp \
	highfive/H5DataType.hpp \
	highfive/H5File.hpp \
	highfive/H5FileDriver.hpp \
	highfive/H5Object.hpp \
	highfive/H5Selection.hpp \
	highfive/H5Attribute.hpp \
	highfive/H5DataSpace.hpp \
	highfive/H5Exception.hpp \
	highfive/H5Group.hpp \
	highfive/H5PropertyList.hpp \
	highfive/H5Utility.hpp \
	highfive/bits/H5Annotate_traits.hpp \
	highfive/bits/H5DataSet_misc.hpp \
	highfive/bits/H5FileDriver_misc.hpp \
	highfive/bits/H5Node_traits.hpp \
	highfive/bits/H5Selection_misc.hpp \
	highfive/bits/H5Annotate_traits_misc.hpp \
	highfive/bits/H5Dataspace_misc.hpp \
	highfive/bits/H5File_misc.hpp \
	highfive/bits/H5Node_traits_misc.hpp \
	highfive/bits/H5Slice_traits.hpp \
	highfive/bits/H5Attribute_misc.hpp \
	highfive/bits/H5DataType_misc.hpp \
	highfive/bits/H5Group_misc.hpp \
	highfive/bits/H5Object_misc.hpp \
	highfive/bits/H5Slice_traits_misc.hpp \
	highfive/bits/H5Converter_misc.hpp \
	highfive/bits/H5Exception_misc.hpp \
	highfive/bits/H5Iterables_misc.hpp \
	highfive/bits/H5PropertyList_misc.hpp \
	highfive/bits/H5Utils.hpp

libSherpaHDF5v1Output_la_LIBADD = @CONDITIONAL_HDF5LIBS@

libSherpaHDF5v1Output_la_CXXFLAGS = $(AM_CXXFLAGS)

libSherpaHDF5v1Output_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_HDF5INCS@

libSherpaHDF5Output_la_LIBADD = @CONDITIONAL_HDF5LIBS@

libSherpaHDF5Output_la_CXXFLAGS = $(AM_CXXFLAGS)

libSherpaHDF5Output_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_HDF5INCS@


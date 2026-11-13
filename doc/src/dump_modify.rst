.. index:: dump\_modify

dump\_modify command
====================

**Syntax:**


.. parsed-literal::

   dump_modify dump-ID keyword values ...

* dump-ID = ID of dump to modify
* one or more keyword/value pairs may be appended
* these keywords apply to various dump styles
* keyword = *append* or *binary* or *buffer* or *every* or *fileper* or *first* or *flush* or *format* or *nfile* or *pad* or *region* or *thresh*
  
  .. parsed-literal::
  
       append arg = yes or no
       binary arg = yes or no
       buffer arg = yes or no
       every arg = N
         N = dump every this many timesteps
         N can be a variable (see below)
       fileper arg = Np
         Np = write one file for every this many processors
       first arg = yes or no
       flush arg = yes or no
       format args = line string, int string, float string, M string, or none
         string = C-style format string
         M = integer from 1 to N, where N = # of per-atom quantities being output
       nfile arg = Nf
         Nf = write this many files, one from each of Nf processors
       pad arg = Nchar = # of characters to convert timestep to
       region arg = region-ID or "none"
       thresh args = attribute operation value
         attribute = same attributes (x,fy,etotal,sxx,etc) used by dump custom style
         operation = "<" or "<=" or ">" or ">=" or "==" or "!="
         value = numeric value to compare to
         these 3 args can be replaced by the word "none" to turn off thresholding

* these keywords apply only to the *image* and *movie* :doc:`styles <dump_image>`
* keyword = *axestrans* or *backcolor* or *backcolor2* or *bitrate* or *boxcolor* or *boxtrans* or *cmap* or *color* or *framerate* or *gamma* or *gcolor* or *glinecolor* or *glinetrans* or *gridgroup* or *gtrans* or *lights* or *metal* or *metalfinish* or *pcolor* or *pdiam* or *ptrans* or *scolor* or *slinecolor* or *slinetrans* or *specular* or *ssaosamples* or *strans* or *subboxcolor* or *subboxtrans* or *surfgroup*
  
  .. parsed-literal::
  
       axestrans arg = opacity
         opacity = opacity of the XYZ axes, from 0.0 (invisible) to 1.0 (opaque)
       backcolor arg = color
         color = name of color for background
       backcolor2 arg = color
         color = name of second color for background gradient, or none
       bitrate arg = rate
         rate = target bitrate for movie in kbps
       boxcolor arg = color
         color = name of color for box lines
       boxtrans arg = opacity
         opacity = opacity of the simulation box lines, from 0.0 to 1.0
       cmap args = mode lo hi style delta N entry1 entry2 ... entryN
         mode = particle or grid or surf or gridx or gridy or gridz
         lo = number or min = lower bound of range of color map
         hi = number or max = upper bound of range of color map
         style = 2 letters = "c" or "d" or "s" plus "a" or "f"
           "c" for continuous
           "d" for discrete
           "s" for sequential
           "a" for absolute
           "f" for fractional
         delta = binsize (only used for style "s", otherwise ignored)
           binsize = range is divided into bins of this width
         N = # of subsequent entries
         entry = value color (for continuous style)
           value = number or min or max = single value within range
           color = name of color used for that value
         entry = lo hi color (for discrete style)
           lo/hi = number or min or max = lower/upper bound of subset of range
           color = name of color used for that subset of values
         entry = color (for sequential style)
           color = name of color used for a bin of values
       color args = name R G B
         name = name of color
         R,G,B = red/green/blue numeric values from 0.0 to 1.0
       framerate arg = fps
         fps = frames per second for movie
       gamma arg = gvalue
         gvalue = gamma adjustment applied to rendered objects (from 0.1 to 10.0, 1.0 = no change)
       gcolor args = proc color
         proc = proc ID or range of IDs (see below)
         color = name of color or color1/color2/...
       glinecolor arg = color
         color = name of color for grid cell outlines
       glinetrans arg = opacity
         opacity = opacity of the grid cell outlines, from 0.0 to 1.0
       gridgroup arg = group-ID
         group-ID = name of a user-defined grid group, see the :doc:`group <group>` command
       gtrans arg = opacity
         opacity = opacity of the grid cells and cutting planes, from 0.0 to 1.0
       lights args = ambient key fill back
         ambient = intensity of ambient light from 0.0 to 1.0
         key = intensity of key light from 0.0 to 1.0
         fill = intensity of fill light from 0.0 to 1.0
         back = intensity of back light from 0.0 to 1.0
       metal arg = mfactor
         mfactor = how metallic the rendered objects appear, from 0.0 (paint) to 1.0 (bare metal)
       metalfinish arg = style
         style = satin or polished or mirror = surface finish of metallic objects
       pcolor args = type color
         type = particle type or range of types or proc ID or range of IDs (see below)
         color = name of color or color1/color2/...
       pdiam args = type diam
         type = particle type or range of types (see below)
         diam = diameter of particles of that type (distance units)
       ptrans args = type opacity
         type = particle type or range of types (see below)
         opacity = opacity of particles of that type, from 0.0 to 1.0
       scolor args = proc color
         proc = proc ID or range of IDs (see below)
         color = name of color for surf one option
       slinecolor arg = color
         color = name of color for surface element outlines
       slinetrans arg = opacity
         opacity = opacity of the surface element outlines, from 0.0 to 1.0
       specular arg = style
         style = none or wide or narrow or tight = specular highlights off or their width
       ssaosamples arg = nsamples
         nsamples = number of SSAO sampling directions per pixel (from 4 to 64)
       strans arg = opacity
         opacity = opacity of the surface elements, from 0.0 to 1.0
       subboxcolor arg = color
         color = name of color for processor sub-box lines
       subboxtrans arg = opacity
         opacity = opacity of the processor sub-box lines, from 0.0 to 1.0
       surfgroup arg = group-ID
         group-ID = name of a user-defined surf group, see the :doc:`group <group>` command



**Examples:**


.. parsed-literal::

   dump_modify 1 format line "%d %d %20.15g %g %g"
   dump_modify 1 format float %20.15g
   dump_modify myDump thresh x < 0.0 thresh vx >= 3.0
   dump_modify 1 every 1000
   dump_modify 1 every v_myVar
   dump_modify 1 cmap particle min max cf 0.0 3 min green 0.5 yellow max blue boxcolor red

**Description:**

Modify the parameters of a previously defined dump command.  Not all
parameters are relevant to all dump styles.





These keywords apply to all dump styles unless otherwise noted.  The
descriptions give details.


----------


The *append* keyword applies to all dump styles except *image* and
*movie*\ .  It also applies only to text output files, not to binary or
gzipped files.  If specified as *yes*\ , then dump snapshots are
appended to the end of an existing dump file.  If specified as *no*\ ,
then a new dump file will be created which will overwrite an existing
file with the same name.  This keyword can only take effect if the
dump\_modify command is used after the :doc:`dump <dump>` command, but
before the first command that causes dump snapshots to be output,
e.g. a :doc:`run <run>` command.  Once the dump file has been opened,
this keyword has no further effect.


----------


The *binary* keyword applies only to the :doc:`dump particle/vtk, dump grid/vtk, and dump surf/vtk <dump_vtk>` styles of the VTK package.  If specified as *yes*\ ,
the VTK data is written in compact binary form; if *no* (the default), it is
written as ASCII which is human-readable but typically requires larger file
sizes than binary. This is a per-dump setting because the VTK styles select
their file format from the filename extension, which would otherwise conflict
with SPARTA's usual ".bin" binary-file convention.


----------


The *buffer* keyword applies only all dump styles except *image* and
*movie*\ .  It also applies only to text output files, not to binary or
gzipped files.  If specified as *yes*\ , which is the default, then each
processor writes its output into an internal text buffer, which is
then sent to the processor(s) which perform file writes, and written
by those processors(s) as one large chunk of text.  If specified as
*no*\ , each processor sends its per-atom data in binary format to the
processor(s) which perform file writes, and those processor(s) format
and write it line by line into the output file.

The buffering mode is typically faster since each processor does the
relatively expensive task of formatting the output for its own atoms.
However it requires about twice the memory (per processor) for the
extra buffering.


----------


The *every* keyword changes the dump frequency originally specified by
the :doc:`dump <dump>` command to a new value.  The every keyword can be
specified in one of two ways.  It can be a numeric value in which case
it must be > 0.  Or it can be an :doc:`equal-style variable <variable>`,
which should be specified as v\_name, where name is the variable name.
In this case, the variable is evaluated at the beginning of a run to
determine the next timestep at which a dump snapshot will be written
out.  On that timestep, the variable will be evaluated again to
determine the next timestep, etc.  Thus the variable should return
timestep values.  See the stagger() and logfreq() math functions for
:doc:`equal-style variables <variable>`, as examples of useful functions
to use in this context.  Other similar math functions could easily be
added as options for :doc:`equal-style variables <variable>`.  When
using the variable option with the *every* keyword, you also need to
use the *first* option if you want an initial snapshot written to the
dump file.

For example, the following commands will
write snapshots at timesteps 0,10,20,30,100,200,300,1000,2000,etc:


.. parsed-literal::

   variable	        s equal logfreq(10,3,10)
   dump		1 particle all 100 tmp.dump id type x y z
   dump_modify	1 every v_s first yes


----------


The *fileper* keyword is documented below with the *nfile* keyword.


----------


The *first* keyword determines whether a dump snapshot is written on
the very first timestep after the dump command is invoked.  This will
always occur if the current timestep is a multiple of N, the frequency
specified in the :doc:`dump <dump>` command, including timestep 0.  But
if this is not the case, a dump snapshot will only be written if the
setting of this keyword is *yes*\ .  If it is *no*\ , which is the
default, then it will not be written.


----------


The *flush* keyword applies to all dump styles except *image* and
*movie*\ .  It also applies only when the styles are used to write
multiple successive snapshots to the same file.  It determines whether
a flush operation is invoked after a dump snapshot is written to the
dump file.  A flush insures the output in that file is current (no
buffering by the OS), even if SPARTA halts before the simulation
completes.


----------


The *format* keyword can be used to change the default numeric format
output by the text-based dump styles: *particle*\ , *grid*\ , *surf*\ ,
*tally*\ .

All the specified format strings are C-style formats, e.g. as used by
the C/C++ printf() command.  The *line* keyword takes a single
argument which is the format string for an entire line of output with
N fields for each particle, grid cell, or suraface elememt, which you
must enclose in quotes if it is more than one field.  The *int* and
*float* keywords take a single format argument and are applied to all
integer or floating-point quantities output.  The setting for *M
string* also takes a single format argument which is used for the Mth
value output in each line, e.g. the 5th column is output in high
precision for "format 5 %20.15g".

The *format* keyword can be used multiple times.  The precedence is
that for each value in a line of output, the *M* format (if specified)
is used, else the *int* or *float* setting (if specified) is used,
else the *line* setting (if specified) for that value is used, else
the default setting is used.  A setting of *none* clears all previous
settings, reverting all values to their default format.

NOTE: Grid cell IDs are stored internally as 4-byte or 8-byte signed
integers, depending on how SPARTA was compiled.  When specifying the
*format int* option you can use a "%d"-style format identifier in the
format string and SPARTA will convert this to the corresponding 8-byte
form it it is needed when outputting those values.  However, when
specifying the *line* option or *format M string* option for those
values, you should specify a format string appropriate for an 8-byte
signed integer, e.g. one with "%ld", if SPARTA was compiled with the
-DSPARTA\_BIGBIG option for 8-byte IDs.


----------


The *nfile* or *fileper* keywords apply to all dump styles except
*image* and *movie*\ .  They can be used in conjunction with the "%"
wildcard character in the specified dump file name.  As explained on
the :doc:`dump <dump>` command doc page, the "%" character causes the
dump file to be written in pieces, one piece for each of P processors.
By default P = the number of processors the simulation is running on.
The *nfile* or *fileper* keyword can be used to set P to a smaller
value, which can be more efficient when running on a large number of
processors.

The *nfile* keyword sets P to the specified Nf value.  For example, if
Nf = 4, and the simulation is running on 100 processors, 4 files will
be written, by processors 0,25,50,75.  Each will collect information
from itself and the next 24 processors and write it to a dump file.

For the *fileper* keyword, the specified value of Np means write one
file for every Np processors.  For example, if Np = 4, every 4th
processor (0,4,8,12,etc) will collect information from itself and the
next 3 processors and write it to a dump file.


----------


The *pad* keyword only applies when the dump filename is specified
with a wildcard "\*" character which becomes the timestep.  If *pad* is
0, which is the default, the timestep is converted into a string of
unpadded length, e.g. 100 or 12000 or 2000000.  When *pad* is
specified with *Nchar* > 0, the string is padded with leading zeroes
so they are all the same length = *Nchar*\ .  For example, pad 7 would
yield 0000100, 0012000, 2000000.  This can be useful so that
post-processing programs can easily read the files in ascending
timestep order.


----------


The *region* keyword only applies to the dump *particle* and *image*
styles.  If specified, only particles in the region will be written to
the dump file or included in the image.  Only one region can be
applied as a filter (the last one specified).  See the
:doc:`region <region>` command for more details.  Note that a region can
be defined as the "inside" or "outside" of a geometric shape, and it
can be the "union" or "intersection" of a series of simpler regions.


----------


The *thresh* keyword only applies to the dump *particle* and *image*
styles.  Multiple thresholds can be specified.  Specifying "none"
turns off all threshold criteria.  If thresholds are specified, only
particles whose attributes meet all the threshold criteria are written
to the dump file or included in the image.  The possible attributes
that can be tested for are the same as those that can be specified in
the :doc:`dump particle <dump>` command.  Note that different attributes
can be output by the dump particle command than are used as threshold
criteria by the dump\_modify command.  E.g. you can output the
coordinates of particles whose velocity components are above some
threshold.





These keywords apply only to the :doc:`dump image <dump_image>` and
:doc:`dump movie <dump_image>` styles.  Any keyword that affects an
image, also affects a movie, since the movie is simply a collection of
images.  Some of the keywords only affect the :doc:`dump movie <dump_image>` style.  The descriptions give details.


----------


The *axestrans*\ , *boxtrans*\ , *glinetrans*\ , *gtrans*\ , *ptrans*\ ,
*slinetrans*\ , *strans* and *subboxtrans* keywords can be used with the
:doc:`dump image <dump_image>` command to render an object as partially
transparent.  Each sets the opacity of one kind of drawn object: the
XYZ axes, the simulation box lines, the grid cell outlines, the grid
cells and cutting planes, the particles, the surface element outlines,
the surface elements, and the processor sub-box lines, respectively.
The opacity value must be between 0.0 (invisible) and 1.0 (fully
opaque); the default for all of them is 1.0.  An opacity of 0.0 renders
exactly as if the object were not drawn at all.

The most common use is to make a surface transparent with *strans* so
that the particles and grid cells inside it become visible, or to thin
out a dense particle field with *ptrans* so that a surface inside it can
be seen.

The *ptrans* keyword sets the opacity per particle species.  The
specified *type* should be an integer from 1 to Nspecies, and a wildcard
asterisk can be used in place of or in conjunction with the *type*
argument to specify a range of species, in the same manner as described
for the *pdiam* keyword.  Note that a per-species setting only takes
effect if the particle color or diameter setting of the :doc:`dump image <dump_image>` command is *type*\ , since only in that case is the
species of each particle known to the dump.  A setting for the full
range of species, e.g. *ptrans \* 0.5*, is applied to all particles
regardless of the color and diameter settings.

Transparency is rendered with the so-called screen-door method: only a
subset of the pixels of an object is drawn, chosen from a repeating
16x16 ordered dither pattern.  No blending of colors is performed.  This
has three consequences worth knowing about.

First, the effect looks best when the *fsaa* keyword of the :doc:`dump image <dump_image>` command is also enabled, since the image is then
rendered at twice the size and averaged down, which blends the drawn and
undrawn pixels together.  Recommended opacity values are 0.25, 0.5, and
0.75.

Second, the dither pattern is anchored to the screen and not to the
object.  Two objects with the *same* opacity that overlap on the screen
therefore leave out the same pixels, so the more distant of the two is
not visible through the nearer one.  Objects only show through each
other if their opacity values differ.

Third, the pixels that are left out are indistinguishable from the
background in the depth buffer, which the *ssao* keyword and the
post-processing of the *outline*\ , *depthcue* and *defocus* keywords of
the :doc:`dump image <dump_image>` command all read.  See the discussion
of those keywords on the :doc:`dump image <dump_image>` doc page.

IMPORTANT NOTE: Transparency should not be combined with the *outline*
keyword of the :doc:`dump image <dump_image>` command.  Outlines are drawn
where the distance to the viewer jumps, and every pixel left out by the
dither pattern looks like such a jump, so a transparent object is
outlined pixel by pixel and is filled almost entirely with the outline
color.


----------


The *backcolor* keyword can be used with the :doc:`dump image <dump_image>` command to set the background color of the
images.  The color name can be any of the 140 pre-defined colors (see
below) or a color name defined by the dump\_modify color option.


----------


The *backcolor2* keyword can be used with the :doc:`dump image <dump_image>` command to set a second background color.  If
set to a color name, the background of the images becomes a vertical
gradient from the *backcolor* color at the bottom of the image to the
*backcolor2* color at the top.  If set to *none*\ , the gradient is
turned off and the uniform *backcolor* color is used.  The color name
can be any of the 140 pre-defined colors (see below) or a color name
defined by the dump\_modify color option.


----------


The *bitrate* keyword can be used with the :doc:`dump movie <dump_image>` command to define the size of the resulting
movie file and its quality via setting how many kbits per second are
to be used for the movie file. Higher bitrates require less
compression and will result in higher quality movies.  The quality is
also determined by the compression format and encoder.  The default
setting is 2000 kbit/s, which will result in average quality with
older compression formats.

IMPORTANT NOTE: Not all movie file formats supported by dump movie
allow the bitrate to be set.  If not, the setting is silently ignored.


----------


The *boxcolor* keyword can be used with the :doc:`dump image <dump_image>` command to set the color of the simulation box
drawn around the particles in each image.  See the "dump image box"
command for how to specify that a box be drawn.  The color name can be
any of the 140 pre-defined colors (see below) or a color name defined
by the dump\_modify color option.


----------


The *cmap* keyword can be used with the :doc:`dump image <dump_image>`
command to define a color map that is used to draw "objects" which can
be particles, grid cells, or surface elements.  The mode setting must
be *particle* or *grid* or *surf* or *gridx* or *gridy* or *gridz* which
correspond to the same keywords in the :doc:`dump image <dump_image>`
command.

Color maps are used to assign a specific RGB (red/green/blue) color
value to an individual object when it is drawn, based on the object's
attribute, which is a numeric value, e.g. the x-component of velocity
for a particle, if the particle-attribute "vx" was specified in the
:doc:`dump image <dump_image>` command.

The basic idea of a color map is that the attribute will be within a
range of values, and that range is associated with a a series of
colors (e.g. red, blue, green).  A specific value (vx = -3.2) can then
mapped to the series of colors (e.g. halfway between red and blue),
and a specific color is determined via an interpolation procedure.

There are many possible options for the color map, enabled by the
*cmap* keyword.  Here are the details.

The *lo* and *hi* settings determine the range of values allowed for
the attribute.  If numeric values are used for *lo* and/or *hi*\ , then
values that are lower/higher than that value are set to the value.
I.e. the range is static.  If *lo* is specified as *min* or *hi* as
*max* then the range is dynamic, and the lower and/or upper bound will
be calculated each time an image is drawn, based on the set of objects
being visualized.

The *style* setting is two letters, such as "ca".  The first letter is
either "c" for continuous, "d" for discrete, or "s" for sequential.
The second letter is either "a" for absolute, or "f" for fractional.

A continuous color map is one in which the color changes continuously
from value to value within the range.  A discrete color map is one in
which discrete colors are assigned to sub-ranges of values within the
range.  A sequential color map is one in which discrete colors are
assigned to a sequence of sub-ranges of values covering the entire
range.

An absolute color map is one in which the values to which colors are
assigned are specified explicitly as values within the range.  A
fractional color map is one in which the values to which colors are
assigned are specified as a fractional portion of the range.  For
example if the range is from -10.0 to 10.0, and the color red is to be
assigned to objects with a value of 5.0, then for an absolute color
map the number 5.0 would be used.  But for a fractional map, the
number 0.75 would be used since 5.0 is 3/4 of the way from -10.0 to
10.0.

The *delta* setting is only specified if the style is sequential.  It
specifies the bin size to use within the range for assigning
consecutive colors to.  For example, if the range is from -10.0 to
10.0 and a *delta* of 1.0 is used, then 20 colors will be assigned to
the range.  The first will be from -10.0 <= color1 < -9.0, then 2nd
from -9.0 <= color2 < -8.0, etc.

The *N* setting is how many entries follow.  The format of the entries
depends on whether the color map style is continuous, discrete or
sequential.  In all cases the *color* setting can be any of the 140
pre-defined colors (see below) or a color name defined by the
dump\_modify color option.

For continuous color maps, each entry has a *value* and a *color*\ .
The *value* is either a number within the range of values or *min* or
*max*\ .  The *value* of the first entry must be *min* and the *value*
of the last entry must be *max*\ .  Any entries in between must have
increasing values.  Note that numeric values can be specified either
as absolute numbers or as fractions (0.0 to 1.0) of the range,
depending on the "a" or "f" in the style setting for the color map.

Here is how the entries are used to determine the color of an
individual object, given the value X of its attribute.  X will fall
between 2 of the entry values.  The color of the object is linearly
interpolated (in each of the RGB values) between the 2 colors
associated with those entries.  For example, if X = -5.0 and the 2
surrounding entries are "red" at -10.0 and "blue" at 0.0, then the
object's color will be halfway between "red" and "blue", which happens
to be "purple".

For discrete color maps, each entry has a *lo* and *hi* value and a
*color*\ .  The *lo* and *hi* settings are either numbers within the
range of values or *lo* can be *min* or *hi* can be *max*\ .  The *lo*
and *hi* settings of the last entry must be *min* and *max*\ .  Other
entries can have any *lo* and *hi* values and the sub-ranges of
different values can overlap.  Note that numeric *lo* and *hi* values
can be specified either as absolute numbers or as fractions (0.0 to
1.0) of the range, depending on the "a" or "f" in the style setting
for the color map.

Here is how the entries are used to determine the color of an
individual object, given the value X of its attribute.  The entries
are scanned from first to last.  The first time that *lo* <= X <=
*hi*\ , X is assigned the color associated with that entry.  You can
think of the last entry as assigning a default color (since it will
always be matched by X), and the earlier entries as colors that
override the default.  Also note that no interpolation of a color RGB
is done.  All objects will be drawn with one of the colors in the list
of entries.

For sequential color maps, each entry has only a *color*\ .  Here is how
the entries are used to determine the color of an individual object,
given the value X of its attribute.  The range is partitioned into N
bins of width *binsize*\ .  Thus X will fall in a specific bin from 1 to
N, say the Mth bin.  If it falls on a boundary between 2 bins, it is
considered to be in the higher of the 2 bins.  Each bin is assigned a
color from the E entries.  If E < N, then the colors are repeated.
For example if 2 entries with colors red and green are specified, then
the odd numbered bins will be red and the even bins green.  The color
of the object is the color of its bin.  Note that the sequential
color map is really a shorthand way of defining a discrete color map
without having to specify where all the bin boundaries are.


----------


The *color* keyword can be used with the :doc:`dump image <dump_image>`
command to define a new color name, in addition to the 140-predefined
colors (see below), and associates 3 red/green/blue RGB values with
that color name.  The color name can then be used with any other
dump\_modify keyword that takes a color name as a value.  The RGB
values should each be floating point values between 0.0 and 1.0
inclusive.

When a color name is converted to RGB values, the user-defined color
names are searched first, then the 140 pre-defined color names.  This
means you can also use the *color* keyword to overwrite one of the
pre-defined color names with new RBG values.


----------


The *framerate* keyword can be used with the :doc:`dump movie <dump_image>` command to define the duration of the resulting
movie file.  Movie files written by the dump *movie* command have a
default frame rate of 24 frames per second and the images generated
will be converted at that rate.  Thus a sequence of 1000 dump images
will result in a movie of about 42 seconds.  To make a movie run
longer you can either generate images more frequently or lower the
frame rate.  To speed a movie up, you can do the inverse.  Using a
frame rate higher than 24 is not recommended, as it will result in
simply dropping the rendered images. It is more efficient to dump
images less frequently.


----------


The *gamma* keyword can be used with the :doc:`dump image <dump_image>` command to adjust the gamma value of the
rendered objects: the summed up light contributions of each pixel are
raised to the power of 1/\ *gvalue* before they are converted to the
8-bit color values of the image file, similar to the gamma adjustment
of image manipulation programs.  A value larger than 1.0 lightens the
image, most strongly in the darker regions, and thus can bring out
shading detail on the dimly lit side of objects; a value smaller than
1.0 darkens the image and increases the contrast.  The default
rendering is already tuned to look right on typical displays, and no
gamma information is stored in the image files, so this setting
fine-tunes the tonal balance rather than applying a required display
correction.  The adjustment applies only to rendered objects; the
background colors are used exactly as specified.


----------


The *gcolor* keyword can be used one or more times with the :doc:`dump image <dump_image>` command, only when its grid color setting is
*proc*\ , to set the color that grid cells will be drawn in the image.

The *proc* setting should be an integer from 1 to Nprocs = the number
of processors.  A wildcard asterisk can be used in place of or in
conjunction with the *proc* argument to specify a range of processor
IDs.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
number of processors, then an asterisk with no numeric values means
all procs from 1 to N.  A leading asterisk means all procs from 1 to n
(inclusive).  A trailing asterisk means all procs from n to N
(inclusive).  A middle asterisk means all procs from m to n
(inclusive).  Note that for this command, processor IDs range from 1
to Nprocs inclusive, instead of the more customary 0 to Nprocs-1.

The specified *color* can be a single color which is any of the 140
pre-defined colors (see below) or a color name defined by the
dump\_modify color option.  Or it can be two or more colors separated
by a "/" character, e.g. red/green/blue.  In the former case, that
color is assigned to all the specified processors.  In the latter
case, the list of colors are assigned in a round-robin fashion to each
of the specified processors.


----------


The *glinecolor* keyword can be used with the :doc:`dump image <dump_image>` command to set the color of the grid cell
outlines drawn around the grid cells in each image.  See the "dump
image gline" command for how to specify that cell outlines be drawn.
The color name can be any of the 140 pre-defined colors (see below) or
a color name defined by the dump\_modify color option.


----------


The *gridgroup* keyword can be used with the :doc:`dump image <dump_image>` command to only draw a subset of the grid cells
in the simulation.  By default all the grid cells are rendered.  The
group-ID argument can be any valid grid group name, as defined by the
:doc:`group grid <group>` command.


----------


The *lights* keyword can be used with the :doc:`dump image <dump_image>` command to adjust the lighting of the rendered
scene.  The scene is illuminated by an ambient light, a key light, a
fill light, and a back light, whose grayscale intensities are set by
the four arguments, each a value from 0.0 (off) to 1.0 (brightest).
The default lighting is ambient 0.0, key 0.9, fill 0.45, back 0.9.


----------


The *metal* keyword can be used with the :doc:`dump image <dump_image>`
command to make objects look like they are made of metal rather than
of colored plastic.  Painted surfaces scatter light in all directions
and reflect it without changing its color, which is what the rendering
does by default; polished metal instead reflects light back directly
and colors it in the process.  An *mfactor* value of 0.0 keeps the
default appearance, 1.0 renders bare metal, and values in between
blend the two.

Since metal shows its surroundings rather than a color of its own, the
rendering approximates them with a bright sky above and a dark ground
below.  This is why metallic objects have a bright upper and a dark
lower half, and why they need a dark background to look convincing.
The color assigned to an object tints these reflections, so gold
objects show a golden sheen and aluminum ones a neutral gray sheen.
The pre-defined color *silver* is a good match for aluminum.  Colors
chosen for paint are often too saturated for metal: for gold, for
example, defining a color with the *color* keyword using the values
1.0 0.766 0.336 looks more like the metal than the pre-defined color
*gold* does.

IMPORTANT NOTE: These settings imitate the appearance of metal with a
few extra operations per pixel.  They do not compute how light
actually travels through the scene.  For images where the appearance
of the material matters, exporting the particle positions and
rendering them with a ray tracing program will always give better
results than any combination of these settings.


----------


The *metalfinish* keyword can be used with the :doc:`dump image <dump_image>` command to select the surface finish used when
the *metal* keyword is enabled.  A *satin* finish has a broad soft
sheen and resembles brushed metal such as aluminum.  A *polished*
finish concentrates the sheen into a narrower streak.  A *mirror*
finish reflects the surroundings the way a curved mirror does, which
makes spheres look like polished ball bearings, with a darker lower
half than the other two settings.  This keyword has no effect unless
the *metal* keyword is set to a value larger than 0.0.


----------


The *pcolor* keyword can be used one or more times with the :doc:`dump image <dump_image>` command, only when its particle color setting is
*type* or *procs*\ , to set the color that particles will be drawn in
the image.

If the particle color setting is *type*\ , then the specified *type* for the
*pcolor* keyword should be an integer from 1 to Ntypes = the number of
particle types.  A wildcard asterisk can be used in place of or in
conjunction with the *type* argument to specify a range of particle
types.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
number of particle types, then an asterisk with no numeric values
means all types from 1 to N.  A leading asterisk means all types from
1 to n (inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).

If the particle color setting is *proc*\ , then the specified *type* for the
*pcolor* keyword should be an integer from 1 to Nprocs = the number of
processors.  A wildcard asterisk can be used in place of or in
conjunction with the *type* argument to specify a range of processor
IDs, just as described above for particle types.  Note that for this
command, processor IDs range from 1 to Nprocs inclusive, instead of
the more customary 0 to Nprocs-1.

The specified *color* can be a single color which is any of the 140
pre-defined colors (see below) or a color name defined by the
dump\_modify color option.  Or it can be two or more colors separated
by a "/" character, e.g. red/green/blue.  In the former case, that
color is assigned to all the specified particle types.  In the latter
case, the list of colors are assigned in a round-robin fashion to each
of the specified particle types.


----------


The *pdiam* keyword can be used with the :doc:`dump image <dump_image>`
command, when its particle diameter setting is *type*\ , to set the size
that particles of each type will be drawn in the image.  The specified
*type* should be an integer from 1 to Ntypes.  As with the *pcolor*
keyword, a wildcard asterisk can be used as part of the *type*
argument to specify a range of particle types.  The specified *diam*
is the size in whatever distance :doc:`units <units>` the input script
is using.


----------


The *scolor* keyword can be used one or more times with the :doc:`dump image <dump_image>` command, only when its surface element color
setting is *one* or *proc*\ , to set the color that surface elements
will be drawn in the image.

When the surf color is *one*\ , the *proc* setting for this command
is ignored.

When the surf color is *proc*\ , the *proc* setting for this command
should be an integer from 1 to Nprocs = the number of processors.  A
wildcard asterisk can be used in place of or in conjunction with the
*proc* argument to specify a range of processor IDs.  This takes the
form "\*" or "\*n" or "n\*" or "m\*n".  If N = the number of processors,
then an asterisk with no numeric values means all procs from 1 to N.
A leading asterisk means all procs from 1 to n (inclusive).  A
trailing asterisk means all procs from n to N (inclusive).  A middle
asterisk means all procs from m to n (inclusive).  Note that for this
command, processor IDs range from 1 to Nprocs inclusive, instead of
the more customary 0 to Nprocs-1.

When the surf color is *one*\ , the specified *color* setting for
this command must be a single color which is any of the 140
pre-defined colors (see below) or a color name defined by the
dump\_modify color option.

When the surf color is *proc*\ , the *color* setting for this command
can be one or more colors separated by a "/" character,
e.g. red/green/blue.  For a single color, that color is assigned to
all the specified processors.  For two or more colors, the list of
colors are assigned in a round-robin fashion to each of the specified
processors.


----------


The *slinecolor* keyword can be used with the :doc:`dump image <dump_image>` command to set the color of the surface element
outlines drawn around the surface elements in each image.  See the
"dump image sline" command for how to specify that surface element
outlines be drawn.  The color name can be any of the 140 pre-defined
colors (see below) or a color name defined by the dump\_modify color
option.


----------


The *specular* keyword can be used with the :doc:`dump image <dump_image>` command to adjust the specular highlights
independently from the *shiny* keyword of the dump image command.  The
*none* setting turns the highlights off entirely, which results in a
rough, matte surface appearance from the remaining diffuse lighting.
The other settings select the width of the highlights: *wide*
highlights are close to the default appearance, *narrow* highlights
are visibly smaller, and *tight* produces small sharp highlights with
a plastic-like appearance.  The *sfactor* value of the *shiny* keyword
scales the brightness of the highlights; without the *specular*
keyword it also sets their width.


----------


The *ssaosamples* keyword can be used with the :doc:`dump image <dump_image>` command to set the number of directions that
the SSAO depth shading enabled by the *ssao* keyword examines around
each pixel.  More directions produce smoother shading; fewer
directions render proportionally faster but make the shading grainier.
Without this setting the number of directions is derived from the
*dfactor* value of the *ssao* keyword and ranges from 8 to 40.
Reducing the number of directions is a simple way to trade some image
quality for faster image output, for example for preview renderings.
The graininess is less visible when the *fsaa* keyword of the dump
image command is also enabled.


----------


The *subboxcolor* keyword can be used with the :doc:`dump image <dump_image>` command to set the color of the processor
sub-box boundaries drawn in each image.  See the "dump image subbox"
command for how to specify that sub-box boundaries be drawn.  The
color name can be any of the 140 pre-defined colors (see below) or a
color name defined by the dump\_modify color option.


----------


The *surfgroup* keyword can be used with the :doc:`dump image <dump_image>` command to only draw a subset of the surface
elements in the simulation.  By default all the surface elements are
rendered.  The group-ID argument can be any valid surf group name, as
defined by the :doc:`group surf <group>` command.


----------


**Restrictions:** none

**Related commands:**

:doc:`dump <dump>`, :doc:`dump image <dump_image>`, :doc:`undump <undump>`

**Default:**

The option defaults are

* append = no
* buffer = yes for all dump styles except *image* and *movie*
* axestrans = 1.0
* backcolor = black
* backcolor2 = none
* boxcolor = yellow
* boxtrans = 1.0
* subboxcolor = yellow
* lights = 0.0 0.9 0.45 0.9
* cmap = mode min max cf 0.0 2 min blue max red, for all modes
* color = 140 color names are pre-defined as listed below
* every = whatever it was set to via the :doc:`dump <dump>` command
* fileper = # of processors
* first = no
* flush = yes
* format = %d and %g for each integer or floating point value
* gamma = 1.0
* gcolor = \* red/green/blue/yellow/aqua/cyan
* glinecolor = white
* glinetrans = 1.0
* gridgroup = all
* gtrans = 1.0
* metal = 0.0
* metalfinish = satin
* nfile = 1
* pad = 0
* pcolor = \* red/green/blue/yellow/aqua/cyan
* pdiam = \* 1.0
* ptrans = \* 1.0
* region = none
* scolor = \* gray
* slinecolor = white
* slinetrans = 1.0
* specular = width derived from the *shiny* keyword of the dump image command
* ssaosamples = number derived from the *dfactor* value of the *ssao* keyword
* strans = 1.0
* subboxtrans = 1.0
* surfgroup = all
* thresh = none


----------


These are the 140 colors that SPARTA pre-defines for use with the
:doc:`dump image <dump_image>` and dump\_modify commands.  Additional
colors can be defined with the dump\_modify color command.  The 3
numbers listed for each name are the RGB (red/green/blue) values.
Divide each value by 255 to get the equivalent 0.0 to 1.0 value.

+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| aliceblue = 240, 248, 255     | antiquewhite = 250, 235, 215         | aqua = 0, 255, 255              | aquamarine = 127, 255, 212     | azure = 240, 255, 255          |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| beige = 245, 245, 220         | bisque = 255, 228, 196               | black = 0, 0, 0                 | blanchedalmond = 255, 255, 205 | blue = 0, 0, 255               |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| blueviolet = 138, 43, 226     | brown = 165, 42, 42                  | burlywood = 222, 184, 135       | cadetblue = 95, 158, 160       | chartreuse = 127, 255, 0       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| chocolate = 210, 105, 30      | coral = 255, 127, 80                 | cornflowerblue = 100, 149, 237  | cornsilk = 255, 248, 220       | crimson = 220, 20, 60          |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| cyan = 0, 255, 255            | darkblue = 0, 0, 139                 | darkcyan = 0, 139, 139          | darkgoldenrod = 184, 134, 11   | darkgray = 169, 169, 169       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| darkgreen = 0, 100, 0         | darkkhaki = 189, 183, 107            | darkmagenta = 139, 0, 139       | darkolivegreen = 85, 107, 47   | darkorange = 255, 140, 0       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| darkorchid = 153, 50, 204     | darkred = 139, 0, 0                  | darksalmon = 233, 150, 122      | darkseagreen = 143, 188, 143   | darkslateblue = 72, 61, 139    |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| darkslategray = 47, 79, 79    | darkturquoise = 0, 206, 209          | darkviolet = 148, 0, 211        | deeppink = 255, 20, 147        | deepskyblue = 0, 191, 255      |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| dimgray = 105, 105, 105       | dodgerblue = 30, 144, 255            | firebrick = 178, 34, 34         | floralwhite = 255, 250, 240    | forestgreen = 34, 139, 34      |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| fuchsia = 255, 0, 255         | gainsboro = 220, 220, 220            | ghostwhite = 248, 248, 255      | gold = 255, 215, 0             | goldenrod = 218, 165, 32       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| gray = 128, 128, 128          | green = 0, 128, 0                    | greenyellow = 173, 255, 47      | honeydew = 240, 255, 240       | hotpink = 255, 105, 180        |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| indianred = 205, 92, 92       | indigo = 75, 0, 130                  | ivory = 255, 240, 240           | khaki = 240, 230, 140          | lavender = 230, 230, 250       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| lavenderblush = 255, 240, 245 | lawngreen = 124, 252, 0              | lemonchiffon = 255, 250, 205    | lightblue = 173, 216, 230      | lightcoral = 240, 128, 128     |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| lightcyan = 224, 255, 255     | lightgoldenrodyellow = 250, 250, 210 | lightgreen = 144, 238, 144      | lightgrey = 211, 211, 211      | lightpink = 255, 182, 193      |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| lightsalmon = 255, 160, 122   | lightseagreen = 32, 178, 170         | lightskyblue = 135, 206, 250    | lightslategray = 119, 136, 153 | lightsteelblue = 176, 196, 222 |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| lightyellow = 255, 255, 224   | lime = 0, 255, 0                     | limegreen = 50, 205, 50         | linen = 250, 240, 230          | magenta = 255, 0, 255          |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| maroon = 128, 0, 0            | mediumaquamarine = 102, 205, 170     | mediumblue = 0, 0, 205          | mediumorchid = 186, 85, 211    | mediumpurple = 147, 112, 219   |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| mediumseagreen = 60, 179, 113 | mediumslateblue = 123, 104, 238      | mediumspringgreen = 0, 250, 154 | mediumturquoise = 72, 209, 204 | mediumvioletred = 199, 21, 133 |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| midnightblue = 25, 25, 112    | mintcream = 245, 255, 250            | mistyrose = 255, 228, 225       | moccasin = 255, 228, 181       | navajowhite = 255, 222, 173    |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| navy = 0, 0, 128              | oldlace = 253, 245, 230              | olive = 128, 128, 0             | olivedrab = 107, 142, 35       | orange = 255, 165, 0           |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| orangered = 255, 69, 0        | orchid = 218, 112, 214               | palegoldenrod = 238, 232, 170   | palegreen = 152, 251, 152      | paleturquoise = 175, 238, 238  |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| palevioletred = 219, 112, 147 | papayawhip = 255, 239, 213           | peachpuff = 255, 239, 213       | peru = 205, 133, 63            | pink = 255, 192, 203           |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| plum = 221, 160, 221          | powderblue = 176, 224, 230           | purple = 128, 0, 128            | red = 255, 0, 0                | rosybrown = 188, 143, 143      |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| royalblue = 65, 105, 225      | saddlebrown = 139, 69, 19            | salmon = 250, 128, 114          | sandybrown = 244, 164, 96      | seagreen = 46, 139, 87         |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| seashell = 255, 245, 238      | sienna = 160, 82, 45                 | silver = 192, 192, 192          | skyblue = 135, 206, 235        | slateblue = 106, 90, 205       |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| slategray = 112, 128, 144     | snow = 255, 250, 250                 | springgreen = 0, 255, 127       | steelblue = 70, 130, 180       | tan = 210, 180, 140            |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| teal = 0, 128, 128            | thistle = 216, 191, 216              | tomato = 253, 99, 71            | turquoise = 64, 224, 208       | violet = 238, 130, 238         |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+
| wheat = 245, 222, 179         | white = 255, 255, 255                | whitesmoke = 245, 245, 245      | yellow = 255, 255, 0           | yellowgreen = 154, 205, 50     |
+-------------------------------+--------------------------------------+---------------------------------+--------------------------------+--------------------------------+


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html

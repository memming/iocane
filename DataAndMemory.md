Internal data types and memory issues

# Introduction #

To be able to process large amount of spike train data, we need to optimize the resources. To avoid repeated calculations, once a reusable quantity is computed it is stored.

# Spike Count #
The number of spikes are in uint32 (0 to 4,294,967,295).
Variable name **N** is reserved for spike counts.

# Spike Timing #
Double precision floating point (64 bit = 8 bytes, 16 decimal digits) is used.
Default values of iocane project assumes the timing is in **seconds**.

# Spike Train #
Since we assume double precision timing and no sampling rate, we cannot use sparse arrays. The alternative is to use the cell array, but cell array takes 60 bytes (or 120 bytes in 64 bit machines). Thus, we use array of all spike times with several indexing arrays to indicate channel, trial, etc. (currently we use cell array)
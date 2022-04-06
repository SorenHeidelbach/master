
from mapped_signal_files import *

nat = HDF5Reader("/shared-nfs/SH/data/zymoHMW/NAT/megalodon/signal_mappings.hdf5")

nat2 = nat._load_reads_batch("Batch_0")

read1 = nat2.get(b"001b4a9e-1766-40b3-8af8-20bf51e5f9d9")


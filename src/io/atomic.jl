# Atomic file replacement shared by the savers.

# Replace `path` atomically: `write!(tmp)` creates and fills a unique temp path in the
# destination's directory (same filesystem, so the file gets the usual umask-derived
# mode), then rename(2) swaps it in. A crash mid-write leaves the existing file
# untouched. rename(2) fails on a directory destination rather than deleting it, and
# replaces a symlink/hardlink as a directory entry instead of writing through it. The
# isdir pre-check only gives a clearer message; rename's own failure is the guarantee.
# The temp is removed on any failure, including InterruptException.
function _replace_atomically(write!::Function, path::AbstractString)
    isdir(path) && throw(ArgumentError("$path is a directory; refusing to replace it"))
    tmp = tempname(dirname(abspath(path)); cleanup = false)
    try
        write!(tmp)
        Base.Filesystem.rename(tmp, path)
    catch
        rm(tmp; force = true)
        rethrow()
    end
    return path
end

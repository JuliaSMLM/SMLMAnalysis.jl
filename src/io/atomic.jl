# Atomic file replacement shared by the savers.

# Replace `path` atomically: `mktempdir` creates an exclusively-created temp
# directory in the destination's directory (same filesystem, so the rename
# below is atomic). Its own mode is not inherited by the file written inside
# it, which keeps the usual umask-derived mode. Writing to a fixed filename
# inside that directory, rather than to a `tempname`-reserved path, means no
# other user sharing the destination directory can pre-place a symlink there
# ahead of us — on POSIX, where the directory is created 0700. On Windows the
# directory inherits the parent's ACL instead, so this protection is POSIX-only
# there (a documented limitation, not a guarantee). `write!(tmp)` fills it, then rename(2) swaps it
# in. A crash mid-write leaves the existing file untouched. rename(2) fails on
# a directory destination rather than deleting it, and replaces a
# symlink/hardlink as a directory entry instead of writing through it. The
# isdir pre-check only gives a clearer message; rename's own failure is the
# guarantee. The temp directory (empty on success, holding the partial write
# on failure) is removed on any exit, including InterruptException.
function _replace_atomically(write!::Function, path::AbstractString)
    isdir(path) && throw(ArgumentError("$path is a directory; refusing to replace it"))
    tmpdir = mktempdir(dirname(abspath(path)); cleanup = false)
    tmp = joinpath(tmpdir, basename(path))
    try
        write!(tmp)
        Base.Filesystem.rename(tmp, path)
    finally
        rm(tmpdir; recursive = true, force = true)
    end
    return path
end

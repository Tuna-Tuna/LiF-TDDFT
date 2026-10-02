function test_native_grids()
% Run from MATLAB/Octave after adding tests/ to the search path.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'matlab','analysis'));
axes_in = {0:3, 0:2, 0:1};
a = reshape(1:24,[2,3,4]);
b = a+2;
[ac,bc,volume] = compare_on_native_grids(a,b,axes_in,axes_in);
assert(isequal(ac,permute(a,[2,3,1])));
assert(isequal(bc,permute(b,[2,3,1])));
assert(volume==1);
bad_axes = axes_in; bad_axes{1}=bad_axes{1}+0.1;
expect_error(@() compare_on_native_grids(a,b,axes_in,bad_axes),'nativegrid:matching');
bad=b; bad(1)=NaN;
expect_error(@() compare_on_native_grids(a,bad,axes_in,axes_in),'nativegrid:missing');
bad_axes=axes_in;bad_axes{1}=[0,1,2,4];
expect_error(@() compare_on_native_grids(a,b,bad_axes,bad_axes),'nativegrid:uniform');
expect_error(@() compare_on_native_grids(a,b(:,:,1:3),axes_in,axes_in),'nativegrid:shape');
disp('Native-grid comparison: 5 behavior checks passed.');
end

function expect_error(fn,expected)
try
    fn();
catch err
    assert(strcmp(err.identifier,expected),err.message);
    return;
end
error('Expected failure was not raised: %s',expected);
end

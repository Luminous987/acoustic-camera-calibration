% extract the offset of the poses and the landmarks

function [poses, landmarks] = get_poses_landmarks(g)

poses = [];
landmarks = [];

for n = 1:length(g.idLookup)

  dim = g.idLookup(n).dimension;
  offset = g.idLookup(n).offset;
  if (offset >= 3*g.M)
    poses = [poses; offset];
  elseif (dim < 3*g.M)
    landmarks = [landmarks; offset];
  end
end

end

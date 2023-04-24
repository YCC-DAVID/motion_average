from pathlib import Path

qinpose_path = Path(r"M:\MoAve2\data_for_moAve_2\ver1\undistorted_Img\Ori_QIN.qin")
qinposev2_path = qinpose_path.with_suffix('.qinv2')

lines = qinpose_path.read_text().splitlines()
v2lines = []
num_cam = int(lines.pop(0))
v2lines.append(str(num_cam))

intrinsics = lines.pop(0).split(' ')
print(intrinsics)

for i in range(num_cam):
  extrinsics = lines.pop(0)
  tokens = extrinsics.split(' ')
  tokens = tokens[0:1] + intrinsics + tokens[1:]
  v2lines.append(' '.join(tokens))

v2text = '\r'.join(v2lines)

qinposev2_path.write_text(v2text)

print('Done!')


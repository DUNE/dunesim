import sys

def read_file(filename):
  with open(filename, 'r') as f:
    lines = [l.split('/')[-1].split() for l in f.readlines()]
  lines = [
    [int(l[0].split('_')[-1].strip('.root')), ' '.join(l)]
    for l in lines
  ]
  lines.sort()
  return lines


def write_lines(filename, lines):
  with open(filename, 'w') as f:
    lines = [': ' + l[1]+'\n' for l in lines]
    f.writelines(lines)


if __name__ == '__main__':
  lines = read_file(sys.argv[1])
  write_lines(sys.argv[2], lines)

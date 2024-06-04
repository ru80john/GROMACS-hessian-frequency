with open('../STRUCTURE/lig.pdb', 'r') as f:
    lines = f.readlines()
    atomlist = {}
    for n, i in enumerate(lines):
        if i == '\n':
            continue
        elif i.split()[0] == 'HETATM':
            if i.split()[2][-1].isdigit():
                atomlist[i.split()[1]] = i.split()[2][:-1]
            else:
                atomlist[i.split()[1]] = i.split()[2]

print('[ bondtypes ]')

bond_const = {}
angle_const = {}

with open('Python_Modified_Scaled_Seminario.sb', 'r') as f:
    lines = f.readlines()
    for index, i in enumerate(lines[1:]):
        index += 1
        if i == '\n':
            continue
        elif '*' in i:
            break
        else:
            value = i.split()
            force_const = float(value[1]) * 100
            bond_length = float(value[2]) * 0.1
            bond_const[value[0]] = f'{bond_length:>12.4f}  {force_const:>.3f}'
    for i in lines[index + 1:]:
        if i == '\n':
            continue
        else:
            value = i.split()
            force_const = float(value[1]) * 4.184
            angle = float(value[2])
            angle_const[value[0]] = f'{angle:>11.3f}{force_const:>11.3f}'

with open('Modified_Seminario_Bonds', 'r') as f:
    lines = f.readlines()
    for i in lines:
        word = i.split()
        print(f'  {atomlist[word[-2]]:<3s}   {atomlist[word[-1]]:<3s}   1' + bond_const[word[0]])

print()
print('[ angletypes ]')
with open('Modified_Seminario_Angle', 'r') as f:
    lines = f.readlines()
    for i in lines:
        word = i.split()
        print(f'  {atomlist[word[-3]]:<3s}   {atomlist[word[-2]]:<3s}   {atomlist[word[-1]]:<3s}    1' + angle_const[word[0]])

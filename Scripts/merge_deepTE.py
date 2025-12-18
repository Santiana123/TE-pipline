import argparse
import re

# 定义替换字典
replace_dict = {
    'ClassII_DNA_Harbinger_MITE': 'DNA/PIF-Harbinger',
    'ClassII_DNA_CACTA_MITE': 'DNA/DTC',
    'ClassII_DNA_TcMar_MITE': 'DNA/DTT',
    'ClassII_DNA_Mutator_MITE': 'DNA/DTM',
    'ClassII_DNA_hAT_MITE': 'hAT',
    'ClassII_MITE': 'DNA',
    'ClassI_LTR_Gypsy': 'LTR/Gypsy',
    'ClassI_LTR_Copia': 'LTR/Copia',
    'ClassI_nLTR_SINE_tRNA': 'SINE',
    'ClassI_nLTR_LINE_L1': 'LINE/L1',
    'unknown': 'Unknown',
    'ClassIII_Helitron': 'RC/Helitron',
    'ClassII': 'DNA',
    'ClassII_DNA_CACTA_Unknown': 'DNA/DTC',
    'ClassII_DNA_CACTA_nMITE': 'DNA/DTC',
    'ClassII_DNA_Harbinger_Unknown': 'DNA/PIF-Harbinger',
    'ClassII_DNA_Harbinger_nMITE': 'DNA/PIF-Harbinger',
    'ClassII_DNA_Mutator_Unknown': 'DNA/DTM',
    'ClassII_DNA_Mutator_nMITE': 'DNA/DTM',
    'ClassII_DNA_P_nMITE': 'DNA',
    'ClassII_DNA_TcMar_Unknown': 'DNA/DTT',
    'ClassII_DNA_TcMar_nMITE': 'DNA/DTT',
    'ClassII_DNA_hAT_Unknown': 'hAT',
    'ClassII_DNA_hAT_nMITE': 'hAT',
    'ClassII_nMITE': 'DNA',
    'ClassI_LTR': 'LTR',
    'ClassI_nLTR': 'Retroposon',
    'ClassI_nLTR_DIRS': 'Retroposon',
    'ClassI_nLTR_PLE': 'Retroposon',
    'ClassI_nLTR_SINE': 'SINE',
    'ClassII_DNA_CACTA_unknown': 'DNA/DTC',
    'ClassII_DNA_Harbinger_unknown': 'DNA/PIF-Harbinger',
    'ClassII_DNA_Mutator_unknown': 'DNA/DTM',
    'ClassII_DNA_TcMar_unknown': 'DNA/DTT',
    'ClassII_DNA_hAT_unknown': 'hAT',
    'ClassI_nLTR_SINE_7SL': 'SINE',
    'ClassI': 'LTR/unknown'
}

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        skip = False
        for line in infile:
            if '#Unknown' in line:
                skip = True  # 如果找到包含#Unknown的行，开始跳过写入
            elif line.strip() == "":  # 如果遇到空行，继续跳过
                continue
            elif not skip:
                # 执行替换操作，确保完全匹配
                for old, new in replace_dict.items():
                    # 使用正则表达式确保完整匹配替换
                    line = re.sub(r'\b' + re.escape(old) + r'\b', new, line)
                outfile.write(line)  # 如果不跳过，则写入
            if line.startswith(">") and skip:  # 如果下一行是新的序列标签，停止跳过
                skip = False
    print(f"处理完成，结果已保存至 {output_file}")

# 主函数，使用argparse来处理命令行参数
def main():
    parser = argparse.ArgumentParser(description='处理DeepTE序列，跳过#Unknown部分并进行替换')
    parser.add_argument('-i', '--input', required=True, help='输入文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')

    args = parser.parse_args()

    process_file(args.input, args.output)

if __name__ == "__main__":
    main()

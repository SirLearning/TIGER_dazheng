#!/bin/bash

# FastCall2 GPU异构优化使用示例
# 本脚本展示如何使用GPU加速的scan步骤

echo "=== FastCall2 GPU异构优化示例 ==="

# 设置基本路径
JAVA_HOME="/usr/lib/jvm/java-11-openjdk"
FASTCALL2_JAR="target/debugFastCall2-1.0-SNAPSHOT.jar"
SAMTOOLS_PATH="/usr/bin/samtools"

# 输入文件路径
REFERENCE_GENOME="data/reference/genome.fa"
TAXA_BAM_MAP="data/input/taxaBamMap.txt"
VARIATION_LIBRARY="data/vLib/variations.vlib"
OUTPUT_DIR="output/gpu_optimized"

# GPU优化配置
GPU_ENABLED="true"          # 启用GPU加速
MIN_GPU_DEPTH="50"          # 使用GPU的最小测序深度
BATCH_SIZE="1000"           # GPU批处理大小
THREADS="32"                # 线程数

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

echo "开始GPU优化的基因型扫描..."
echo "GPU加速: $GPU_ENABLED"
echo "最小GPU深度: $MIN_GPU_DEPTH"
echo "批处理大小: $BATCH_SIZE"
echo "线程数: $THREADS"
echo ""

# 执行FastCall2 scan步骤（GPU优化版本）
$JAVA_HOME/bin/java -Xmx64g -Djava.library.path=/usr/local/cuda/lib64 \
    -jar "$FASTCALL2_JAR" \
    -app FastCall2 \
    -mod scan \
    -a "$REFERENCE_GENOME" \
    -b "$TAXA_BAM_MAP" \
    -c "$VARIATION_LIBRARY" \
    -d "1:1,50000000" \
    -e 0 \
    -f 30 \
    -g 20 \
    -h 0.05 \
    -i "$SAMTOOLS_PATH" \
    -j "$THREADS" \
    -k "$OUTPUT_DIR" \
    -gpu "$GPU_ENABLED"

echo ""
echo "GPU优化的基因型扫描完成！"
echo "结果保存在: $OUTPUT_DIR"

# 显示性能统计
if [ -f "$OUTPUT_DIR/performance_report.txt" ]; then
    echo ""
    echo "=== 性能报告 ==="
    cat "$OUTPUT_DIR/performance_report.txt"
fi

echo ""
echo "=== 优化效果对比 ==="
echo "建议运行以下命令进行性能对比："
echo "1. CPU版本: -gpu false"
echo "2. GPU版本: -gpu true"
echo "3. 比较两者的运行时间和资源使用情况"

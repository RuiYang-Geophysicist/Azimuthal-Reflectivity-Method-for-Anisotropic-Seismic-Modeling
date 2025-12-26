#!/bin/bash
# 替换 GitHub 仓库的完整指令脚本
# 目标仓库: https://github.com/RuiYang-Geophysicist/Azimuthal-Reflectivity-Method-for-Anisotropic-Seismic-Modeling.git

set -e  # 遇到错误立即退出

echo "=========================================="
echo "  替换 GitHub 仓库内容"
echo "=========================================="
echo ""

# 1. 检查当前状态
echo "步骤 1: 检查 Git 状态..."
git status

echo ""
read -p "确认要继续吗？(y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "操作已取消"
    exit 1
fi

# 2. 添加所有更改（包括删除的文件）
echo ""
echo "步骤 2: 添加所有更改到暂存区..."
git add -A

# 3. 提交更改
echo ""
echo "步骤 3: 提交更改..."
git commit -m "Clean up: Remove MATLAB files and unnecessary build artifacts

- Removed all MATLAB .m source files (replaced by C++/Python)
- Removed test/comparison scripts
- Removed build artifacts (.so files)
- Removed development notes (Memo.md)
- Removed test data files (.mat, .png)

Project now focuses on C++/Python implementation only."

# 4. 显示提交信息
echo ""
echo "步骤 4: 查看提交历史..."
git log --oneline -5

# 5. 推送到远程仓库
echo ""
echo "步骤 5: 推送到远程仓库..."
echo "警告: 这将覆盖远程仓库的内容！"
read -p "确认要强制推送吗？(y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "推送已取消。你可以稍后手动运行: git push origin main"
    exit 0
fi

# 强制推送（替换远程仓库）
git push origin main --force

echo ""
echo "=========================================="
echo "  ✓ 仓库替换完成！"
echo "=========================================="
echo ""
echo "仓库地址: https://github.com/RuiYang-Geophysicist/Azimuthal-Reflectivity-Method-for-Anisotropic-Seismic-Modeling.git"
echo ""


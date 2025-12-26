# 如何分享私有仓库给测试人员

## 方法一：添加协作者（Collaborators）- 推荐

### 步骤：

1. **打开仓库设置**
   - 访问：https://github.com/RuiYang-Geophysicist/Azimuthal-Reflectivity-Method-for-Anisotropic-Seismic-Modeling
   - 点击仓库页面右上角的 **Settings**（设置）

2. **进入协作者管理**
   - 在左侧菜单中找到 **Collaborators**（协作者）
   - 或者直接访问：
     ```
     https://github.com/RuiYang-Geophysicist/Azimuthal-Reflectivity-Method-for-Anisotropic-Seismic-Modeling/settings/access
     ```

3. **添加协作者**
   - 点击 **Add people**（添加人员）按钮
   - 输入测试人员的 GitHub 用户名或邮箱
   - 选择权限级别：
     - **Read**（读取）：只能查看和克隆代码
     - **Triage**（分类）：可以管理 issues 和 pull requests
     - **Write**（写入）：可以推送代码和创建分支
     - **Maintain**（维护）：可以管理仓库设置（推荐给测试人员：**Read** 或 **Write**）

4. **发送邀请**
   - 点击 **Add [username] to this repository**
   - 测试人员会收到邮件通知
   - 他们需要接受邀请才能访问仓库

### 权限说明：

- **Read（读取）**：适合只需要查看和测试代码的用户
- **Write（写入）**：适合需要提交 bug 报告或修复的用户
- **Maintain（维护）**：适合需要管理仓库的协作者

---

## 方法二：使用 GitHub CLI（命令行）

如果你安装了 GitHub CLI，可以使用命令：

```bash
# 添加协作者（读取权限）
gh repo add-collaborator RuiYang-Geophysicist/Azimuthal-Reflectivity-Method-for-Anisotropic-Seismic-Modeling USERNAME --permission read

# 添加协作者（写入权限）
gh repo add-collaborator RuiYang-Geophysicist/Azimuthal-Reflectivity-Method-for-Anisotropic-Seismic-Modeling USERNAME --permission write
```

---

## 方法三：创建临时访问令牌（Token）

如果测试人员不需要长期访问，可以创建临时访问令牌：

1. **生成 Personal Access Token**
   - GitHub Settings → Developer settings → Personal access tokens → Tokens (classic)
   - 点击 **Generate new token**
   - 选择权限：`repo`（完整仓库访问）
   - 设置过期时间（如 30 天）
   - 复制生成的 token

2. **分享给测试人员**
   - 他们可以使用 token 克隆仓库：
     ```bash
     git clone https://TOKEN@github.com/RuiYang-Geophysicist/Azimuthal-Reflectivity-Method-for-Anisotropic-Seismic-Modeling.git
     ```

⚠️ **注意**：Token 一旦生成就无法再次查看，请妥善保管。

---

## 推荐方案

对于测试人员，建议使用 **方法一（添加协作者）**，权限设置为 **Read** 或 **Write**：
- ✅ 简单易用
- ✅ 可以随时撤销访问
- ✅ 测试人员可以正常使用 Git 工作流
- ✅ 可以查看谁在访问仓库

---

## 撤销访问

如果需要移除某个测试人员的访问权限：
1. 进入 Settings → Collaborators
2. 找到要移除的用户
3. 点击用户名旁边的 **X** 或选择 **Remove**（移除）


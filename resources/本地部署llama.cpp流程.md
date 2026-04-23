>系统配置：![[Pasted image 20260424021451.png]]

---
### 第一步：准备 Fedora 编译环境
Fedora 默认工具链很全，先安装构建所需的开发库：
```bash
sudo dnf groupinstall "Development Tools"
sudo dnf install cmake gcc-c++ git wget
```

### 第二步：安装 Intel OneAPI (驱动 B580 的核心)
由于 Intel Arc 显卡在 Linux 上使用 SYCL 后端，你需要 Intel 的编译环境。
1. **添加 Intel 仓库**：
```bash
sudo tee /etc/yum.repos.d/oneAPI.repo <<EOF
[oneAPI]
name=Intel® oneAPI repository
baseurl=https://yum.repos.intel.com/oneapi
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
EOF
```
2. **安装基础运行依赖**：
```bash
sudo dnf install intel-oneapi-compiler-dpcpp-cpp intel-oneapi-mkl intel-oneapi-mkl-devel
# Fedora 将 Level-Zero 驱动整合在了 intel-compute-runtime 及其子包中
sudo dnf install intel-compute-runtime
```
### 第三步：下载并编译 llama.cpp
创建项目目录并进入，然后执行：
```bash
git clone https://github.com/ggerganov/llama.cpp
cd llama.cpp

# 1. 激活 Intel 编译环境 (必做)
source /opt/intel/oneapi/setvars.sh

# 2. 编译
rm -rf build out # 清理旧缓存，如果有的话
export CC=icx
export CXX=icpx

# 开启 SYCL 后端与 F16 计算加速
cmake -B build -DGGML_SYCL=ON -DGGML_SYCL_F16=ON
cmake --build build --config Release -j$(nproc)

# 提取文件到 out 目录 
cmake --install build --prefix ./out
```
4. 赋予权限能够调用独显
```bash
# 赋予用户组权限，需要重启生效
sudo usermod -aG video,render $USER
```
4. 测试是否能够识别独显
```bash
# 1. 设置路径与驱动协议 
export LD_LIBRARY_PATH=$PWD/out/lib64:$LD_LIBRARY_PATH 
export ONEAPI_DEVICE_SELECTOR=opencl:gpu 
# 2. 运行识别工具 
./out/bin/llama-ls-sycl-device
```
正常识别会输出的内容能看到：
`[opencl:gpu:0]| Intel Arc B580 Graphics`

### 第四步：下载模型
1. 先创建目录用来放模型，我的模型放在了llama.cpp/Models中，进入这个目录，然后创建python虚拟环境
```bash
cd ~/llama.cpp/Models
# 1. 创建名为 'ai_env' 的虚拟环境
python3 -m venv ai_env

# 2. 激活环境 (激活后 shell 前缀会多出一个 (ai_env))
source ai_env/bin/activate

# 3. 安装huggingface_hub和hf_transfer用于下载模型
pip install -U pip  # 先升级 pip 自身
pip install huggingface_hub hf_transfer
```
2. 下载（注意下载是关闭代理）
```bash
# 设置加速环境变量
export HF_ENDPOINT=https://hf-mirror.com
export HF_HUB_ENABLE_HF_TRANSFER=1

# 开启多线程并行下载
./ai_env/bin/hf download unsloth/Qwen3.5-35B-A3B-GGUF \
  --include "Qwen3.5-35B-A3B-Q4_K_M.gguf" \
  --local-dir .
```

### 第五步：启动
```bash
# 1. 设置环境变量（如果刚才的终端没关，这步可以跳过）
export LD_LIBRARY_PATH=$PWD/out/lib64:$LD_LIBRARY_PATH
export ONEAPI_DEVICE_SELECTOR=opencl:gpu
export ZES_ENABLE_SYSMAN=1

# 2. 设定层数和上下文大小，实测30层会爆显存
./out/bin/llama-cli -m ./Models/Qwen3.5-35B-A3B-Q4_K_M.gguf \
  --n-gpu-layers 20 \
  --ctx-size 4096 \
  --threads 12 \
  -cnv \
  --chat-template-kwargs '{"enable_thinking":true}' \
  -p "你是什么模型"
```

### 第六步：设置快捷启动脚本
```bash
cat <<EOF > chat.sh
#!/bin/bash
source /opt/intel/oneapi/setvars.sh
export LD_LIBRARY_PATH=\$PWD/out/lib64:\$LD_LIBRARY_PATH
export ONEAPI_DEVICE_SELECTOR=opencl:gpu
export ZES_ENABLE_SYSMAN=1

./out/bin/llama-cli -m ./Models/Qwen3.5-35B-A3B-Q4_K_M.gguf \\
  --n-gpu-layers 25 \\
  --ctx-size 4096 \\
  --threads 12 \\
  -cnv \\
  --chat-template-kwargs '{"enable_thinking":true}' \\
  -p "你是一个有用的智能助手，回答简洁高效，处理问题专业高级。"
EOF

chmod +x chat.sh
```
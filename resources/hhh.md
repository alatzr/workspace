```markdown
# 六合彩投注俚语解析系统 · 完整项目文档 (AI Agent 实现版)

> **版本**: v1.0  
> **定位**: 封闭领域 · 高准确率 · 确定性语言编译器  
> **受众**: AI 编程 Agent（大模型）  
> **目标**: 基于本文档，Agent 可独立完成系统各模块的设计、编码、测试与集成

---

## 目录

1. [项目总览与设计哲学](#1-项目总览与设计哲学)
2. [系统架构](#2-系统架构)
3. [DSL 设计规范](#3-dsl-设计规范)
4. [Layer 1: 词法层 (Tokenizer)](#4-layer-1-词法层-tokenizer)
5. [Layer 2: 语义切分层 (Chunker)](#5-layer-2-语义切分层-chunker)
6. [Layer 3: T5 归一化层](#6-layer-3-t5-归一化层)
7. [Layer 4: DSL 规则解析层 (Parser)](#7-layer-4-dsl-规则解析层-parser)
8. [Layer 5: 校验层 (Validator)](#8-layer-5-校验层-validator)
9. [训练数据工程](#9-训练数据工程)
10. [数据勘探与噪声建模](#10-数据勘探与噪声建模)
11. [评估体系](#11-评估体系)
12. [Pipeline 串联与错误传播](#12-pipeline-串联与错误传播)
13. [项目路线图](#13-项目路线图)
14. [开发者指南与常见陷阱](#14-开发者指南与常见陷阱)
15. [附录](#15-附录)

---

## 1. 项目总览与设计哲学

### 1.1 项目目标

构建一个**高正确率的六合彩投注俚语解析系统**，将真实聊天中的投注文本自动转换为结构化投注数据。

### 1.2 核心设计原则

> **不让 AI 直接负责最终语义。**

系统的核心思想是**分层解耦**：

- **AI (T5) 只负责语言归一化**：风格统一、省略补全、错别字修复、非标准表达转换
- **规则系统负责语义确定**：Chunker 切分、DSL Parser 解析、Validator 校验

这不是一个聊天 AI 项目，而是一个**领域语言编译器 (Domain-Specific Language Compiler)**。

### 1.3 关键约束

- **封闭领域**：仅处理六合彩投注相关俚语，拒绝闲聊
- **100% 确定性**：Parser 和 Validator 不允许任何猜测或概率性输出
- **可追溯**：每层输出必须包含足够元数据，支持回溯到原始文本
- **硬边界**：任何层遇到无法处理的输入，必须报错而非静默丢弃

### 1.4 项目资源

| 资源 | 数量/型号 | 用途 |
|------|-----------|------|
| 真实聊天记录 | 20,000 条 | 数据勘探、噪声建模、盲测 |
| GPU | RTX 3090 20GB | T5 微调训练 |
| 黄金测试集 | 500 条 (人工标注) | 端到端评估 |

---

## 2. 系统架构

### 2.1 流水线总览

```

原始聊天消息
↓
[Layer 1] Tokenizer (词法层)
↓ token 序列
[Layer 2] Chunker (语义切分层)
↓ chunk 列表 (BET / INCOMPLETE / SUMMARY / NOISE)
[Layer 3] T5 Normalizer (归一化层)
↓ DSL 字符串
[Layer 4] DSL Parser (规则解析层)
↓ 结构化投注对象
[Layer 5] Validator (校验层)
↓ 最终结构化 JSON (带校验标记)

```

### 2.2 各层职责矩阵

| 层 | 负责 | 不负责 | 输入 | 输出 |
|----|------|--------|------|------|
| Tokenizer | 字符标准化、token 分类 | 纠错、语义理解 | 原始字符串 | token 序列 |
| Chunker | 金额作用域切分、总额识别 | 字段补全、模式推断 | token 序列 | chunk 列表 |
| T5 Normalizer | 风格统一、省略补全、纠错 | 语义解析、金额校验 | chunk 文本 | DSL 字符串 |
| DSL Parser | DSL → 结构化对象 | 玩法合法性、总额校验 | DSL 字符串 | 投注对象 |
| Validator | 合法性校验、总额反算 | 归一化、切分 | 投注对象列表 | 最终 JSON |

### 2.3 系统非功能需求

- **延迟**：单条消息端到端处理 < 500ms (含 T5 推理)
- **可扩展性**：新增玩法只需扩展词典和 DSL TYPE 表，无需改动架构
- **可调试性**：每层中间结果均可输出，支持逐层回溯
- **高准确率目标**：端到端完全匹配率 > 95% (在黄金测试集上)

---

## 3. DSL 设计规范

### 3.1 格式定义

```

TYPE|TARGETS|MODE|AMOUNT[|SCORE]

```

- 字段分隔符：英文竖线 `|`
- 目标分隔符：英文逗号 `,`
- 禁止任何字段内出现空格
- 第 5 字段 `SCORE` 为可选置信度 (0~1 浮点数)，Parser 和 Validator 忽略

### 3.2 TYPE 编码表

| TYPE | 全称 | 说明 | 最少目标数 | 最多目标数 |
|------|------|------|-----------|-----------|
| `TM` | 特码 | 猜一个特定生肖或号码 | 1 | 无上限 |
| `PT` | 平特 | 不限位置，指定生肖开出即中 | 1 | 无上限 |
| `PTO` | 平特一肖 | 平特的单生肖版 | 1 | 1 |
| `L2` | 二连肖 | 连续两个生肖为一组 | 2 | 无上限 |
| `L3` | 三连肖 | 连续三个生肖为一组 | 3 | 无上限 |
| `BOX` | 包肖 | 包某个/些生肖的所有组合 | 1 | 无上限 |

扩展规则：新增玩法只需在此表增加一行，Parser 的 TYPE 枚举同步更新。

### 3.3 MODE 编码表

| MODE | 含义 | 金额作用对象 | 触发条件 |
|------|------|------------|----------|
| `SINGLE` | 仅一个目标 | 该目标自身 | TARGETS 长度 = 1 |
| `EACH` | 目标列表中每一个独立目标 | 每个目标个体 | TARGETS 长度 ≥ 1 |
| `GROUP` | 给定目标列表即为最终组合 | 整个目标列表的组合 | 连肖类，用户已指定组合 |
| `PERM` | 目标列表中所有符合玩法的排列/组合 | 每一个由列表生成的组合 | 连肖类，用户给出集合需展开 |

#### `GROUP` vs `PERM` 示例

```

输入: 两连肖猴马200          → L2|猴,马|GROUP|200     (用户指定只买"猴-马"这一组)
输入: 两连肖龙马虎200        → L2|龙,马,虎|PERM|200   (生成龙-马、龙-虎、马-虎三组，每组200)

```

### 3.4 TARGETS 规范

- 元素类型：生肖单字 (鼠牛虎兔龙蛇马羊猴鸡狗猪) 或两位数字 (01~49)
- 允许混合：`蛇,09,马` 合法
- 元素分隔：英文逗号 `,`，前后无空格
- 禁止嵌套、禁止分组括号

### 3.5 AMOUNT 规范

- 纯正整数
- 不含货币单位 (元、米、蚊、块等) — 由 T5 归一化层清洗
- 不含货币符号 (¥、$ 等)

### 3.6 DSL 合法性自检规则

以下情况 DSL 视为非法，Parser 必须报错：
- 字段数 ≠ 4 且 ≠ 5
- TYPE 不在预定义集合中
- TARGETS 为空或包含非法字符
- MODE 不在 {SINGLE, EACH, GROUP, PERM} 中
- AMOUNT 不是纯正整数
- TARGETS 长度 = 1 但 MODE ≠ SINGLE

### 3.7 DSL 示例全集

```

基础特码

TM|蛇,狗,鸡|EACH|5
TM|09,33,45|EACH|5
TM|蛇|SINGLE|10

平特

PT|鼠|SINGLE|1500
PT|牛,虎|EACH|200

连肖 GROUP

L2|猴,马|GROUP|300
L3|鼠,牛,虎|GROUP|500

连肖 PERM

L2|龙,马,虎|PERM|200
L3|鼠,牛,虎,兔|PERM|100

包肖

BOX|猪,狗|PERM|50

带置信度

TM|蛇,狗,鸡|EACH|5|0.98

```

---

## 4. Layer 1: 词法层 (Tokenizer)

### 4.1 职责定义

> 将原始字符串转换为语义 token 序列。**不做理解，只做分类标注。**

### 4.2 Token 类型枚举

```python
TOKEN_TYPE = [
    "ANIMAL",      # 生肖: 鼠牛虎兔龙蛇马羊猴鸡狗猪
    "NUMBER",      # 数字: 01~49 或纯数字串
    "PLAY_TYPE",   # 玩法: 特码,平特,平特一肖,两连肖,二连肖,三连肖,包肖
    "MODE",        # 模式词: 各,各数,每个数,包
    "AMOUNT",      # 金额数字 (上下文无关暂不区分，由 Chunker 消歧)
    "SUM_TRIGGER", # 总额触发词: 共,计,合计,总,总额
    "SYMBOL",      # 有意义的分隔符: , . / -
    "UNKNOWN"      # 无法识别但保留的字符
]
```

注意: AMOUNT 类型在当前阶段与 NUMBER 无法区分，Tokenzier 统一标为 NUMBER。AMOUNT 类型保留供未来扩展。

4.3 字符标准化流水线

严格按以下顺序执行：

```
输入字符串
  → Step 1: Emoji 映射
  → Step 2: 全角转半角
  → Step 3: 中文标点转英文标点
  → Step 4: 移除所有空白字符
  → 输出标准化字符串
```

Step 1: Emoji 映射表

```python
EMOJI_MAP = {
    # 生肖 emoji
    "🐭": "鼠", "🐮": "牛", "🐯": "虎", "🐰": "兔",
    "🐲": "龙", "🐍": "蛇", "🐴": "马", "🐏": "羊",
    "🐵": "猴", "🐔": "鸡", "🐶": "狗", "🐷": "猪",
    # 数字 emoji (含变体)
    "0️⃣": "0", "1️⃣": "1", "2️⃣": "2", "3️⃣": "3", "4️⃣": "4",
    "5️⃣": "5", "6️⃣": "6", "7️⃣": "7", "8️⃣": "8", "9️⃣": "9",
    # 不可见变体选择器 (直接删除)
    "️": "", "⃣": ""
}
```

Step 2: 全角转半角映射表

```python
FULL_TO_HALF = {
    # 数字 (必须优先于字母)
    "０": "0", "１": "1", "２": "2", "３": "3", "４": "4",
    "５": "5", "６": "6", "７": "7", "８": "8", "９": "9",
    # 大写字母
    "Ａ": "A", "Ｂ": "B", "Ｃ": "C", "Ｄ": "D", "Ｅ": "E",
    "Ｆ": "F", "Ｇ": "G", "Ｈ": "H", "Ｉ": "I", "Ｊ": "J",
    "Ｋ": "K", "Ｌ": "L", "Ｍ": "M", "Ｎ": "N", "Ｏ": "O",
    "Ｐ": "P", "Ｑ": "Q", "Ｒ": "R", "Ｓ": "S", "Ｔ": "T",
    "Ｕ": "U", "Ｖ": "V", "Ｗ": "W", "Ｘ": "X", "Ｙ": "Y", "Ｚ": "Z",
    # 小写字母
    "ａ": "a", "ｂ": "b", "ｃ": "c", "ｄ": "d", "ｅ": "e",
    "ｆ": "f", "ｇ": "g", "ｈ": "h", "ｉ": "i", "ｊ": "j",
    "ｋ": "k", "ｌ": "l", "ｍ": "m", "ｎ": "n", "ｏ": "o",
    "ｐ": "p", "ｑ": "q", "ｒ": "r", "ｓ": "s", "ｔ": "t",
    "ｕ": "u", "ｖ": "v", "ｗ": "w", "ｘ": "x", "ｙ": "y", "ｚ": "z",
    # 标点符号
    "，": ",", "。": ".", "／": "/", "；": ";", "：": ":",
    "！": "!", "＠": "@", "＃": "#", "＄": "$", "％": "%",
    "＾": "^", "＆": "&", "＊": "*", "（": "(", "）": ")",
    "～": "~", "　": " ", "？": "?", "【": "[", "】": "]",
    "「": "{", "」": "}", "｜": "|", "＼": "\\"
}
```

Step 3: 额外中文标点替换

```python
CN_TO_EN_PUNCT = {
    "、": ",",
    "“": '"', "”": '"',
    "‘": "'", "’": "'",
    "《": "<", "》": ">",
    "…": "..."
}
```

Step 4: 空白字符清除

```python
import re
# 移除所有 Unicode 空白字符
text = re.sub(r'\s+', '', text)
```

4.4 关键词词典 (贪婪匹配)

词典优先级 (按此顺序构建 Trie)

最长匹配优先原则：若 "平特一肖" 命中，则不再匹配 "平特"。

所有关键词列表：

```python
DICTIONARY = {
    "PLAY_TYPE": [
        "平特一肖",   # 必须在 "平特" 前面 (更长)
        "两连肖",     # 别名
        "二连肖",     # 别名
        "三连肖",
        "平特",
        "特码",
        "包肖"
    ],
    "MODE": [
        "每个数",
        "各数",
        "各",
        "包"
    ],
    "SUM_TRIGGER": [
        "合计",
        "总额",
        "总共",
        "共",
        "计",
        "总"
    ],
    "ANIMAL": [
        "鼠", "牛", "虎", "兔", "龙", "蛇",
        "马", "羊", "猴", "鸡", "狗", "猪"
    ]
}
```

Trie 构建伪代码

```python
class TrieNode:
    def __init__(self):
        self.children = {}
        self.token_type = None
        self.is_end = False

class KeywordTrie:
    def __init__(self):
        self.root = TrieNode()
    
    def insert(self, word, token_type):
        """插入关键词，设置 token_type"""
        node = self.root
        for char in word:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.is_end = True
        node.token_type = token_type
    
    def match_longest(self, text, start_pos):
        """从 start_pos 开始做最长匹配，返回 (token_type, length) 或 (None, 0)"""
        node = self.root
        longest_match = None
        longest_len = 0
        for i in range(start_pos, len(text)):
            char = text[i]
            if char not in node.children:
                break
            node = node.children[char]
            if node.is_end:
                longest_match = node.token_type
                longest_len = i - start_pos + 1
        return longest_match, longest_len
```

4.5 Token 分类流程

```python
def tokenize(standardized_text: str) -> List[Tuple[str, str]]:
    """
    输入: 标准化后的文本
    输出: [(token_text, token_type), ...]
    """
    tokens = []
    i = 0
    while i < len(standardized_text):
        char = standardized_text[i]
        
        # 1. 尝试关键词 Trie 匹配 (最长匹配)
        token_type, length = trie.match_longest(standardized_text, i)
        if token_type:
            tokens.append((standardized_text[i:i+length], token_type))
            i += length
            continue
        
        # 2. 数字匹配
        if char.isdigit():
            j = i
            while j < len(standardized_text) and standardized_text[j].isdigit():
                j += 1
            tokens.append((standardized_text[i:j], "NUMBER"))
            i = j
            continue
        
        # 3. 符号匹配
        if char in {',', '.', '/', '-'}:
            tokens.append((char, "SYMBOL"))
            i += 1
            continue
        
        # 4. 未识别字符 (保留)
        tokens.append((char, "UNKNOWN"))
        i += 1
    
    return tokens
```

4.6 Tokenizer 输出示例

```
输入: "🐍狗鸡5米"
标准化: "蛇狗鸡5米"
输出:
  [("蛇", ANIMAL), ("狗", ANIMAL), ("鸡", ANIMAL), ("5", NUMBER), ("米", UNKNOWN)]

输入: "平特一肖猴1500蚊"
标准化: "平特一肖猴1500蚊"
输出:
  [("平特一肖", PLAY_TYPE), ("猴", ANIMAL), ("1500", NUMBER), ("蚊", UNKNOWN)]
```

4.7 Tokenizer 必须遵守的约束

· ❌ 禁止纠错 ("平码" 不能改成 "平特"，只能标 UNKNOWN)
· ❌ 禁止推断语义 ("5" 不能标 AMOUNT，只标 NUMBER)
· ❌ 禁止删除任何 token (UNKNOWN 也必须保留)
· ❌ 禁止重排序
· ❌ 禁止跨 token 合并或拆分

---

5. Layer 2: 语义切分层 (Chunker)

5.1 职责定义

基于金额作用域和玩法边界，将 token 序列切分为独立的投注表达片段。

5.2 设计核心：有限状态机

Chunker 必须用状态转移表驱动，禁止超过两层的 if-else 嵌套。

5.3 Chunk 类型

```python
CHUNK_TYPE = {
    "BET":         "完整投注表达 (含目标和金额)",
    "INCOMPLETE":  "缺失字段 (缺目标或缺金额)，不丢弃，留给下游处理",
    "SUMMARY":     "总额说明行 (含总额触发词)",
    "NOISE":       "纯噪音 (无任何投注相关 token)"
}
```

5.4 状态定义

状态 含义 进入条件
START 初始状态，无待处理数据 系统启动 / chunk 输出后重置
WAIT_TARGET 已有目标列表，等待金额或模式词 START 读到 ANIMAL/NUMBER
WAIT_MODE 已有目标，刚读到模式词，等待金额 WAIT_TARGET 读到 MODE
HAVE_AMOUNT 金额前置，等待目标 START 或 WAIT_TARGET 闭合后读到 AMOUNT 且无目标
IN_TYPE 刚读到玩法词，等待目标列表 START 读到 PLAY_TYPE
SUMMARY_START 读到总额触发词，后续直接累积文本 任意状态读到 SUM_TRIGGER
COMPLETE chunk 已闭合，输出并重置 满足闭合条件
ERROR 非法状态，记录并尝试恢复 非法转移

5.5 状态转移表

Chunker 内部维护以下上下文变量：

· pending_targets: List[str] — 累积的目标 (ANIMAL/NUMBER)
· pending_amount: Optional[int] — 暂存的金额
· current_play_type: Optional[str] — 当前玩法
· current_mode: Optional[str] — 当前模式词
· chunk_type: str — 当前 chunk 类型 (默认 BET)
· output_chunks: List[Dict] — 输出积累

状态转移表 (必须严格实现)：

当前状态 Token 类型 条件 动作 下一状态
START ANIMAL / NUMBER - 加入 pending_targets WAIT_TARGET
START PLAY_TYPE - 记录为 current_play_type IN_TYPE
START NUMBER pending_targets 为空 暂存为 pending_amount；chunk_type=BET HAVE_AMOUNT
START SUM_TRIGGER - chunk_type=SUMMARY；累积文本开始 SUMMARY_START
START MODE - 非法，记录错误 ERROR
START SYMBOL / UNKNOWN - 忽略 START
WAIT_TARGET ANIMAL / NUMBER - 加入 pending_targets WAIT_TARGET
WAIT_TARGET MODE - 记录为 current_mode WAIT_MODE
WAIT_TARGET NUMBER pending_targets 非空 闭合当前 chunk (BET)；清空 pending_targets；暂存此金额为 pending_amount HAVE_AMOUNT
WAIT_TARGET PLAY_TYPE pending_targets 非空 闭合当前 chunk (BET)；重置；记录新 PLAY_TYPE IN_TYPE
WAIT_TARGET PLAY_TYPE pending_targets 为空 替换 current_play_type IN_TYPE
WAIT_TARGET SUM_TRIGGER pending_targets 非空 先闭合当前 chunk (BET)；chunk_type=SUMMARY；开始累积 SUMMARY_START
WAIT_TARGET SUM_TRIGGER pending_targets 为空 chunk_type=SUMMARY；开始累积 SUMMARY_START
WAIT_TARGET SYMBOL - 忽略 (分隔目标用) WAIT_TARGET
WAIT_TARGET UNKNOWN - 加入 pending_targets (保留噪音，留给 T5) WAIT_TARGET
WAIT_MODE NUMBER - 闭合当前 chunk (BET) COMPLETE→START
WAIT_MODE 其他 - 非法，标记当前 chunk 为 INCOMPLETE；输出；重置 ERROR→START
HAVE_AMOUNT ANIMAL / NUMBER - 加入 pending_targets；闭合当前 chunk (BET) COMPLETE→START
HAVE_AMOUNT PLAY_TYPE - 金额前置但无目标：先输出一个 INCOMPLETE chunk (仅有金额)；重置；进入新 chunk IN_TYPE
HAVE_AMOUNT SUM_TRIGGER - 先输出 INCOMPLETE chunk (仅有金额)；切换 SUMMARY SUMMARY_START
HAVE_AMOUNT 其他 - 非法，等待恢复 ERROR
IN_TYPE ANIMAL / NUMBER - 加入 pending_targets WAIT_TARGET
IN_TYPE PLAY_TYPE - 替换 current_play_type (上一个玩法无目标，丢弃) IN_TYPE
IN_TYPE 其他 - 非法 (玩法后必须有目标) ERROR
SUMMARY_START 任意 - 累积文本；直到 token 流结束 SUMMARY_START
ERROR 任意 - 强制闭合当前残 chunk 为 INCOMPLETE；重置；重新处理当前 token START

5.6 Chunk 闭合条件 (优先级从高到低)

1. SUM_TRIGGER：出现总额触发词 → 当前 chunk 闭合，后续全部归为 SUMMARY
2. PLAY_TYPE：出现新玩法词，且 pending_targets 非空 → 当前 chunk 闭合
3. 金额后置正常闭合：WAIT_TARGET 状态 + 读到 NUMBER + pending_targets 非空
4. 金额前置正常闭合：HAVE_AMOUNT 状态 + 读到 ANIMAL/NUMBER
5. 流结束：剩余 pending 数据 → 标记 INCOMPLETE 输出

5.7 输出结构

```python
{
    "chunk_id": int,         # 递增 ID
    "text": str,             # 原始文本片段
    "tokens": [              # token 序列
        {"text": str, "type": str},
        ...
    ],
    "chunk_type": "BET" | "INCOMPLETE" | "SUMMARY" | "NOISE",
    "span": [start, end],    # 原始文本中的字符位置 [闭区间]
    "has_play_type": bool,   # chunk 内是否包含玩法词
    "has_mode": bool,        # chunk 内是否包含模式词
    "has_amount": bool,      # chunk 内是否包含金额
    "has_targets": bool      # chunk 内是否包含目标
}
```

5.8 必须通过的测试用例

# 输入 期望 Chunk 输出 (类型: 文本)
1 猴鼠各数20，猪马龙狗各数10，平特猴800 BET:猴鼠各数20; BET:猪马龙狗各数10; BET:平特猴800
2 23,07,43,44各数20 BET:23,07,43,44各数20
3 蛇狗鸡5 BET:蛇狗鸡5 (无 MODE)
4 20块蛇狗鸡 BET:20块蛇狗鸡 (金额前置)
5 蛇狗平特鸡800 INCOMPLETE:蛇狗; BET:平特鸡800
6 共370 SUMMARY:共370
7 蛇各5狗10鸡各5 BET:蛇各5; BET:狗10; BET:鸡各5
8 平特猴800总计1000 BET:平特猴800; SUMMARY:总计1000
9 蛇,狗,鸡各5 BET:蛇,狗,鸡各5
10 两连肖猴马300 BET:两连肖猴马300

5.9 Chunker 必须遵守的约束

· ❌ 禁止丢弃任何 chunk (INCOMPLETE 也必须输出)
· ❌ 禁止在 chunk 内部做纠错或补全
· ❌ 禁止合并相邻 chunk
· ❌ 禁止猜测缺失的模式词或玩法

---

6. Layer 3: T5 归一化层

6.1 职责定义

将 Chunker 产出的每个 BET chunk 归一化为标准 DSL 字符串。

6.2 模型选择

项目 选择
主力模型 Langboat/mengzi-t5-base (220M 参数)
备选模型 Langboat/mengzi-t5-large (~1B 参数)
对比模型 google/mt5-base (用于词表分析对比)
选择理由 中文优化词表, 生肖/数字 token 完整, 生成更干净

6.3 任务定义

输入: Chunker 输出的 BET chunk 文本 (含噪声、错别字、省略)
输出: 标准 DSL 字符串 (严格符合第 3 节规范)

归一化操作清单:

· 错别字纠正: 平码 → 平特
· 省略补全: 蛇狗鸡5 → TM|蛇,狗,鸡|EACH|5 (补 TYPE 和 MODE)
· 金额后缀清洗: 5米 → 5
· 标点统一: 中文逗号 → 英文逗号
· 空格清除
· 玩法别名统一: 两连肖 → 二连肖 → 但 DSL TYPE 用 L2

6.4 训练数据格式

```json
{
  "source": "蛇狗鸡5米",
  "target": "TM|蛇,狗,鸡|EACH|5"
}
```

6.5 训练超参数推荐

```python
training_args = {
    "model_name": "Langboat/mengzi-t5-base",
    "max_source_length": 128,
    "max_target_length": 64,
    "learning_rate": 5e-5,
    "optimizer": "AdamW",
    "scheduler": "cosine",
    "batch_size": 32,       # 3090 20GB 可轻松支持
    "gradient_accumulation_steps": 1,
    "epochs": 10,
    "warmup_steps": 500,
    "weight_decay": 0.01,
    "dropout": 0.1,
    "fp16": False,          # 3090 显存充足，用全精度
    "early_stopping_patience": 3,
    "eval_steps": 200,
    "save_steps": 500
}
```

6.6 推理参数

```python
generation_args = {
    "max_length": 64,
    "num_beams": 4,         # 束搜索
    "do_sample": False,     # 禁用采样，强制确定性
    "temperature": 0.1,     # 极低温度 (需要时)
    "early_stopping": True,
    "repetition_penalty": 1.2  # 防止重复生成
}
```

6.7 输出后处理

T5 生成后，必须检查：

· 生成的字符串是否满足 DSL 格式
· 若不满足，尝试简单修复 (补齐分隔符)
· 若无法修复，标记为 LOW_CONFIDENCE，输出原始生成内容 + 警告

---

7. Layer 4: DSL 规则解析层 (Parser)

7.1 职责定义

将 DSL 字符串确定性解析为结构化投注对象。100% 确定性，不允许猜测。

7.2 输入输出

```python
# 输入: str
"TM|蛇,狗,鸡|EACH|5"

# 输出: dict
{
    "bet_type": "TM",
    "targets": ["蛇", "狗", "鸡"],
    "mode": "EACH",
    "amount": 5,
    "confidence": None  # 若有第5字段则填入
}
```

7.3 解析规则

```python
def parse_dsl(dsl_string: str) -> dict:
    # Step 1: 按 | 分割
    fields = dsl_string.strip().split('|')
    if len(fields) not in {4, 5}:
        raise ParseError(f"字段数错误: 期望 4 或 5, 实际 {len(fields)}")
    
    bet_type, targets_str, mode, amount_str = fields[:4]
    confidence = float(fields[4]) if len(fields) == 5 else None
    
    # Step 2: TYPE 枚举校验
    VALID_TYPES = {"TM", "PT", "PTO", "L2", "L3", "BOX"}
    if bet_type not in VALID_TYPES:
        raise ParseError(f"未知玩法类型: {bet_type}")
    
    # Step 3: 拆分 TARGETS
    targets = targets_str.split(',')
    if not targets or targets == ['']:
        raise ParseError("目标列表为空")
    
    # Step 4: 验证每个 target 是合法生肖或两位数字
    VALID_ANIMALS = {"鼠","牛","虎","兔","龙","蛇","马","羊","猴","鸡","狗","猪"}
    for t in targets:
        t = t.strip()
        if t in VALID_ANIMALS:
            continue
        if t.isdigit() and 1 <= int(t) <= 49:
            continue
        raise ParseError(f"非法目标: {t}")
    
    # Step 5: MODE 枚举校验
    VALID_MODES = {"SINGLE", "EACH", "GROUP", "PERM"}
    if mode not in VALID_MODES:
        raise ParseError(f"未知模式: {mode}")
    
    # Step 6: AMOUNT 解析
    if not amount_str.isdigit():
        raise ParseError(f"金额非正整数: {amount_str}")
    amount = int(amount_str)
    if amount <= 0:
        raise ParseError(f"金额必须为正: {amount}")
    
    # Step 7: 一致性校验
    if len(targets) == 1 and mode != "SINGLE":
        raise ParseError(f"单目标时 MODE 必须为 SINGLE, 实际: {mode}")
    
    return {
        "bet_type": bet_type,
        "targets": targets,
        "mode": mode,
        "amount": amount,
        "confidence": confidence
    }
```

7.4 PERM 展开算法

当 mode == "PERM" 时，Parser 只保留原始信息，展开由 Validator 完成：

```python
def expand_perm(bet: dict) -> List[dict]:
    """
    输入: 含 PERM 模式的投注对象
    输出: 展开后的投注列表
    """
    from itertools import combinations
    
    bet_type = bet["bet_type"]
    targets = bet["targets"]
    amount = bet["amount"]
    
    if bet_type == "L2":
        combo_size = 2
    elif bet_type == "L3":
        combo_size = 3
    elif bet_type == "BOX":
        combo_size = len(targets)  # 全组合
    else:
        raise ValueError(f"{bet_type} 不支持 PERM 模式")
    
    combos = list(combinations(targets, combo_size))
    
    return [
        {
            "bet_type": bet_type,
            "targets": list(combo),
            "mode": "GROUP",
            "amount": amount,
            "confidence": bet.get("confidence")
        }
        for combo in combos
    ]
```

---

8. Layer 5: 校验层 (Validator)

8.1 职责定义

对解析后的投注对象进行合法性校验，并执行总额反算。

8.2 校验规则列表

规则 ID 规则名称 说明 严重级别
V01 号码范围 NUMBER 必须在 01~49 ERROR
V02 生肖合法 ANIMAL 必须在 12 生肖内 (Parser 已保证) ERROR
V03 金额正整数 金额 > 0 且为整数 (Parser 已保证) ERROR
V04 目标数下限 L2 至少 2 目标, L3 至少 3 目标 ERROR
V05 目标数上限 TM 不加限制; PTO 必须为 1 ERROR
V06 重复投注 同一玩法下完全相同的目标组合出现多次 WARNING
V07 总额校验 所有 BET 金额之和 = SUMMARY 中的总额 ERROR
V08 金额异常大 单注金额 > 100000 WARNING
V09 玩法-模式兼容 单目标只能用 SINGLE；连肖用 GROUP/PERM ERROR
V10 空投注列表 整条消息无有效 BET chunk WARNING

8.3 校验输出结构

```python
{
    "bets": [                    # 所有有效投注 (PERM 已展开)
        {
            "bet_type": "TM",
            "targets": ["蛇"],
            "mode": "SINGLE",
            "amount": 5
        },
        ...
    ],
    "summary": {                 # 总额信息 (若有)
        "text": "共15",
        "claimed_total": 15,
        "calculated_total": 15,
        "match": True
    },
    "errors": [                  # 校验错误列表
        {
            "rule_id": "V06",
            "severity": "WARNING",
            "message": "重复投注: TM 蛇 EACH",
            "detail": "..."
        }
    ],
    "warnings": [...],
    "is_valid": True,            # 无 ERROR 则为 True
    "original_message": "...",   # 原始消息 (追溯)
    "pipeline_metadata": {...}   # 各层耗时、版本等
}
```

8.4 总额反算算法

```python
def validate_total(bets: list, summary_chunks: list) -> dict:
    """
    将所有 BET 的金额总和与 SUMMARY 中声称的总额比对
    """
    calculated = sum(b["amount"] for b in bets)
    
    results = []
    for summary in summary_chunks:
        claimed = extract_total_from_text(summary["text"])
        match = (claimed == calculated)
        results.append({
            "text": summary["text"],
            "claimed_total": claimed,
            "calculated_total": calculated,
            "match": match
        })
    
    return results
```

---

9. 训练数据工程

9.1 数据来源

来源 数量 用途
规则生成器 3,000~5,000 条 T5 训练集
真实数据 (人工标注) 500 条 黄金测试集 (端到端评估)
真实数据 (未标注) 19,500 条 盲测、噪声分析、主动学习

9.2 规则生成器设计

核心思路

```
[模板选择] → [填入生肖/号码] → [应用噪声模型] → 产出 (source, target) 对
```

模板结构定义

```python
TEMPLATES = [
    # (玩法, 目标数范围, 模式, 金额范围, 结构变体)
    {
        "play_type": "TM",
        "targets_count": (1, 6),
        "modes": ["SINGLE", "EACH"],
        "amount_range": (5, 500),
        "variants": [
            "{targets}{mode}{amount}元",           # 标准
            "{amount}元{targets}",                  # 金额前置
            "{amount}元{targets}{mode}",            # 金额前置+模式后置
            "买{targets}{mode}{amount}元",          # 带前缀
            "{targets}{mode}{amount}",              # 无后缀
        ]
    },
    ...
]
```

噪声注入参数

噪声注入参数必须来自数据勘探 (见第 10 节)：

```python
NOISE_PARAMS = {
    "amount_suffix": {
        "enable": True,
        "prob": 0.25,                    # 25% 的样本加金额后缀
        "suffixes": {"米": 0.5, "蚊": 0.2, "块": 0.2, "元": 0.1}  # 概率分布
    },
    "mode_omission": {
        "enable": True,
        "prob": 0.35                     # 35% 省略模式词
    },
    "comma_replacement": {
        "enable": True,
        "prob": 0.40,                    # 40% 用中文逗号或直接连写
        "styles": {"，": 0.4, "": 0.6}   # 中文逗号 vs 连写
    },
    "typo_injection": {
        "enable": True,
        "prob": 0.05,                    # 5% 概率注入错别字
        "typo_map": {
            "平特": "平码",
            "包肖": "包消",
            "连肖": "连消"
        }
    },
    "prefix_injection": {
        "enable": True,
        "prob": 0.15,
        "prefixes": ["买", "下", "投", "跟", "追"]
    }
}
```

---

10. 数据勘探与噪声建模

10.1 目标

从 20,000 条真实聊天记录中自动提取：

· 脏话指纹 (未登录词、错别字、特殊后缀)
· 结构模式分布 (金额位置、标点习惯、省略率)
· 参数化噪声模型的精确概率分布

10.2 勘探步骤

Step 1: 词典掩码分析

```
算法:
  对每条消息:
    用 Tokenizer 词典匹配已知 token
    收集剩余碎片 (未命中词典的连续字符)
  聚合所有碎片 → 按频率排序 → 输出脏话候选清单
```

产出示例:

碎片 频率 归类 处理建议
米 3,240 金额后缀 加入 NOISE_PARAMS
蚊 1,820 金额后缀 (粤语) 加入 NOISE_PARAMS
块 850 金额后缀 加入 NOISE_PARAMS
平码 320 错别字→平特 加入 typo_map
买 2,100 动词前缀 加入前缀列表
下 1,500 动词前缀 加入前缀列表
蚊香 45 不明 (可能闲聊) 标记供人工审核

Step 2: 结构骨架提取

定义结构骨架:

```
{玩法位置: head/tail/none} | {模式词有无: yes/no} | {金额位置: head/tail/middle} | {标点类型: none/comma/cn_comma/mixed} | {前缀有无: yes/no} | {后缀有无: yes/no}
```

对每条消息提取骨架，统计占比：

产出示例:

```
骨架分布 (Top 5):
  none|yes|tail|comma|no|no     → 45%  (如: 蛇,狗,鸡各5)
  head|yes|tail|none|no|no     → 18%  (如: 特码蛇狗鸡各5)
  none|no|tail|none|no|yes     → 12%  (如: 蛇狗鸡5米)
  none|yes|tail|cn_comma|no|no → 10%  (如: 蛇，狗，鸡各5)
  head|yes|head|none|no|no     → 5%   (如: 5元蛇狗鸡各)
```

这些概率直接用于规则生成器的模板权重。

Step 3: 金额与标点专项统计

金额位置统计:

```
尾置: 72%
前置: 18%
中置 (模式词后): 10%
```

标点使用统计:

```
无标点直接连写: 50%
英文逗号: 30%
中文逗号: 15%
混合: 5%
```

Step 4: 碎片合并成噪声参数

将以上统计转换为 NOISE_PARAMS 配置。

10.3 分布比对验证

用以下方法验证生成数据与真实数据分布一致：

```python
def distribution_similarity(real_messages, generated_messages):
    """
    比较两个数据集的:
    - Token 类型分布 (KL 散度)
    - 结构骨架分布 (卡方检验)
    - 消息长度分布 (KS 检验)
    """
    pass
```

---

11. 评估体系

11.1 评估层级

层级 指标 说明 目标
L1: Chunker 边界准确率 chunk 边界完全正确的比例 98%
L2: T5 DSL 字段级准确率 TYPE/TARGETS/MODE/AMOUNT 全部正确 95%
L3: Parser 解析成功率 DSL → 结构化对象无错误 99%
L4: Validator 校验通过率 合法输入通过校验的比例 90%
L5: 端到端 最终 JSON 完全匹配 与人工标注完全一致 92%

11.2 错误分类

```python
ERROR_CATEGORIES = {
    "E_CHUNK_BOUNDARY": "Chunker 切分边界错误",
    "E_T5_MISSING_FIELD": "T5 生成漏字段",
    "E_T5_WRONG_TYPE": "T5 生成 TYPE 错误",
    "E_T5_WRONG_TARGET": "T5 生成 TARGET 错误 (错字/遗漏/多增)",
    "E_T5_WRONG_MODE": "T5 生成 MODE 错误",
    "E_T5_WRONG_AMOUNT": "T5 金额错误 (数字错/错位)",
    "E_T5_FORMAT_ERROR": "T5 生成 DSL 格式不合法",
    "E_PARSE_ERROR": "Parser 解析失败",
    "E_VALIDATION_ERROR": "Validator 校验不通过",
    "E_TOTAL_MISMATCH": "总额不匹配"
}
```

11.3 黄金测试集构建规范

· 从 20,000 条中分层抽样 500 条 (覆盖所有玩法、所有模式、金额前置/后置、多投注、总额行)
· 人工标注每条消息的最终 JSON 输出
· 标注原则：若原文有歧义，选择最合理的解释并备注
· 这 500 条永远不参与训练

---

12. Pipeline 串联与错误传播

12.1 统一数据流对象

```python
class PipelineContext:
    """贯穿整个 Pipeline 的上下文"""
    def __init__(self, original_message: str):
        self.original = original_message
        self.metadata = {
            "pipeline_version": "1.0",
            "timestamp": None,
            "layers": {}  # 每层耗时
        }
        self.tokens = []          # Layer 1 输出
        self.chunks = []          # Layer 2 输出
        self.dsl_strings = []     # Layer 3 输出
        self.bets_raw = []        # Layer 4 输出 (展开前)
        self.bets_expanded = []   # Layer 4 输出 (展开后)
        self.validation_result = None  # Layer 5 输出
        self.errors = []          # 跨层错误收集
```

12.2 错误传播规则

```
Layer 1 Tokenizer 无法分类 → 标记 UNKNOWN，继续传递
Layer 2 Chunker 产出 INCOMPLETE → 传递给 Layer 3，T5 尝试修复
Layer 3 T5 生成非法 DSL → 标记低置信度，传递给 Layer 4
Layer 4 Parser 解析失败 → 错误记录到 errors[]，该 chunk 丢弃
Layer 5 Validator 校验失败 → 错误记录，该 bet 标记为 INVALID
```

12.3 最终输出格式

```python
{
    "status": "success" | "partial" | "error",
    "bets": [...],              # 所有有效投注
    "errors": [...],            # 所有错误
    "warnings": [...],          # 所有警告
    "metadata": {
        "pipeline_version": "1.0",
        "total_time_ms": 123,
        "layers": {
            "tokenizer_ms": 1,
            "chunker_ms": 2,
            "t5_inference_ms": 100,
            "parser_ms": 1,
            "validator_ms": 5
        }
    }
}
```

---

13. 项目路线图

Phase 1: 基础设施 (1-2 天)

· 实现 Tokenizer (tokenizer.py + 词典 + 标准化)
· 实现 Chunker (chunker.py + 状态转移表)
· 实现 Parser (parser.py + PERM 展开)
· 实现 Validator (validator.py + 校验规则)
· 编写所有单元测试

Phase 2: 数据勘探 (1-2 天)

· 对 20,000 条真实数据进行词典掩码分析
· 提取脏话指纹 (碎片频率统计)
· 提取结构骨架分布
· 生成 NOISE_PARAMS 配置
· 生成数据勘探报告

Phase 3: 训练数据生成 (1 天)

· 实现规则生成器
· 挂载噪声模型
· 生成 3,000~5,000 条训练数据
· 分布比对验证 (真实 vs 生成)

Phase 4: T5 微调 (1-2 天)

· 加载 Mengzi-T5-base
· 训练数据格式转换
· 训练 (10 epochs, early stopping)
· 保存最佳模型
· 推理测试

Phase 5: 黄金测试集构建 (1-2 天)

· 从真实数据中分层抽样 500 条
· 人工标注最终 JSON
· 建立评估脚本

Phase 6: 端到端集成 (1 天)

· 实现 Pipeline 串联
· 错误传播机制
· 端到端测试

Phase 7: 迭代优化 (持续)

· 运行黄金测试集评估
· 收集错误案例
· 分析错误原因
· 针对性修复 (补充训练数据/调整规则/修正噪声参数)
· 重复评估循环

---

14. 开发者指南与常见陷阱

14.1 关键设计决策理由

决策 理由
AI 只做归一化，不做语义 生成模型天然有幻觉，投注系统一个字段错误就是致命错误
DSL 作为中间语言 比 JSON 更稳定、易训练、易 debug
状态机驱动 Chunker 比正则或 split 更稳健，可处理金额前置/省略/混用
PERM vs GROUP 消除连肖排列歧义，Parser 不做猜测
硬边界设计 非法输入直接报错，强制上游质量
真实数据不直接训练 无标注；用噪声模型模拟真实分布更可控

14.2 常见实现陷阱

陷阱 表现 避免方法
Tokenizer 最短匹配 "平特一肖" 被切成 "平特" + "一肖" 使用 Trie 确保最长匹配
Chunker 深层 if-else 代码不可维护 使用状态转移表 + 状态模式
训练数据太干净 T5 在真实脏话上泛化差 噪声模型参数必须来自真实数据统计
分词器词表不友好 生肖被切碎，生成质量差 选择 Mengzi-T5，验证 tokenizer 行为
忽视标点 多余空格导致 DSL 解析失败 T5 后处理必须 strip 空格
总额反算遗漏 SUMMARY chunk 被当 BET 处理 Chunker 优先识别 SUM_TRIGGER
Chunker 丢弃残 chunk 有效信息被静默丢弃 INCOMPLETE 必须传递到下游
金额后缀未清洗 Parser 报错 "金额非整数" T5 训练数据中必须包含后缀清洗样本

14.3 调试技巧

· 逐层可视化：每个 Pipeline 阶段都输出中间结果，定位问题所在层
· 错误日志含回溯：每条错误必须包含原始文本 + 层名 + 具体原因
· 小批量金标测试：开发阶段用 10~20 条手动构造的边界样本快速验证

---

15. 附录

15.1 生肖与号码对照 (参考)

```
鼠: 01, 13, 25, 37, 49
牛: 02, 14, 26, 38
虎: 03, 15, 27, 39
兔: 04, 16, 28, 40
龙: 05, 17, 29, 41
蛇: 06, 18, 30, 42
马: 07, 19, 31, 43
羊: 08, 20, 32, 44
猴: 09, 21, 33, 45
鸡: 10, 22, 34, 46
狗: 11, 23, 35, 47
猪: 12, 24, 36, 48
```

15.2 依赖库

```txt
torch>=2.0.0
transformers>=4.30.0
datasets>=2.12.0
scikit-learn>=1.2.0
pandas>=2.0.0
numpy>=1.24.0
wandb>=0.15.0          # 训练可视化 (可选)
tensorboard>=2.13.0    # 训练可视化 (可选)
```

15.3 文件结构建议

```
project_root/
├── data/
│   ├── raw_messages.txt           # 20,000 条原始消息
│   ├── gold_test_set.json        # 500 条人工标注
│   ├── noise_params.json         # 噪声模型参数
│   └── exploration_report.md     # 数据勘探报告
├── src/
│   ├── tokenizer.py              # Layer 1
│   ├── chunker.py                # Layer 2
│   ├── t5_normalizer.py          # Layer 3
│   ├── parser.py                 # Layer 4
│   ├── validator.py              # Layer 5
│   ├── pipeline.py               # 串联
│   ├── data_generator.py         # 规则生成器 + 噪声
│   └── utils.py                  # 工具函数
├── tests/
│   ├── test_tokenizer.py
│   ├── test_chunker.py
│   ├── test_parser.py
│   ├── test_validator.py
│   └── test_pipeline.py
├── models/
│   └── t5_normalizer/            # 微调后的模型
├── configs/
│   └── training_config.json
├── requirements.txt
└── README.md
```

---

文档版本: v1.0
最后更新: 2026-05-08
核心原则: 确定性、可追溯、分层解耦、硬边界

```
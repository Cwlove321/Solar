# 使用说明：预处理后逐阶段运行

操作顺序为：**预处理磁图和轴线 → 生成案例 → PFSS → 背景松弛 → 插绳/加载 → MHD演化**。
每个阶段单独运行，检查完成后再进入下一阶段。

计算环境、依赖库、各阶段编译和作业脚本的准备见[DEPENDENCIES.md](DEPENDENCIES.md)。
以下提交示例采用LSF；使用前应确认当前案例各目录已有编译好的`amrvac`，且`job`
中的环境、队列和核数适合实际机器。案例生成命令只准备源码与数据，不自动编译。

## 1. 预处理磁图和轴线

按[Preprocess/README.md](Preprocess/README.md)选择自己的磁图、范围和轴线，得到：

```text
boundary_product/
  parameters.json
  pfss_boundary.dat
  cutout.dat              # 未生成Brll0000.fits时使用
  Brll0000.fits           # 也可在生成案例时自动生成
axis_product/
  path.dat
  path_summary.json       # 新预处理生成的坐标与参考时刻记录
```

边界需要同一参考时刻的HMI field、inclination、azimuth、disambig分段及全球径向磁图。
轴线可在AIA/HMI图上选取，应与磁图对应同一参考时刻和区域。
`path.dat`为小端float64三元组：`(Stonyhurst经度+180度, 纬度度, 半径/Rsun)`。

## 2. 确认参数，生成案例

```bash
cp case.delivery.example.json science.json
```

按本事件修改`science.json`中的磁绳半径、磁通、螺度、暗条加载、AMR和各阶段时间。
`thermal_profile`只选择04阶段的热物理配置。时间字段以`_code`结尾，换算到秒须乘
启动日志的`unit_time`（当前单位约5972.477 s）。选图和选轴不会自动决定这些科学参数。

```bash
python workflow.py new \
  --boundary /data/my_event/boundary_product \
  --axis /data/my_event/axis_product/path.dat \
  --settings science.json \
  --output /data/runs/my_event

python workflow.py check /data/runs/my_event --config-only
```

程序从实际输入推导角域、角向根网格和轴线点数，同步各阶段参数，并生成对应的
`.par.readme`。输入形状、方向、单位、轴线范围或网格整除不符合要求时会报错。
输出必须是新目录；原始预处理数据和工作流模板不会被修改。

## 3. 按目录逐个运行

以下每段命令**分开执行**。提交上一阶段后，先等待它结束并检查日志和DAT，再运行下一段。
`NNNN`是占位符，必须替换为实际达到该阶段目标时间、且已检查有限性的DAT编号，
不要固定套用0020或0011。若已有`restart.dat`，先核对其来源，不要直接覆盖。

### PFSS初态

```bash
cd /data/runs/my_event/01_pfss
bsub < job
```

若`pfss_stage=01_pfss_global_ct`，改在`01_pfss_global_ct`目录运行；后续链接也改为
该目录输出的`pfss_global_ct_0000.dat`。两个PFSS分支只运行所选的一个。

### 无磁绳背景松弛

```bash
cd /data/runs/my_event/02_relax
ln -s ../01_pfss/data/pfss_0000.dat restart.dat
bsub < job
```

确认初态正常、松弛达到设定时间且背景状态可接受，再选择其输出供插绳使用。

### 插绳及所选质量加载

```bash
cd /data/runs/my_event/03_insert
ln -s ../02_relax/data/relax_NNNN.dat restart.dat
bsub < job
```

该阶段首次从背景DAT启动时执行插绳和所选加载。不要把它自身已插绳的DAT当作背景，
再次按首次插绳设置运行，否则可能重复插绳/加载。插绳后直接保存还是继续加载，
由本案例的质量加载方式和03阶段时间决定。

### 正式MHD演化

```bash
cd /data/runs/my_event/04_fullmhd
ln -s ../03_insert/data/rbsl_insert_load_NNNN.dat restart.dat
bsub < job
```

若配置`fullmhd_source=relax`，这是不插绳的分支：跳过03，04直接链接02的已验收DAT；
生成案例时PFSS也会关闭磁绳Br扣除，并要求不加载暗条。

每次交接至少确认：任务已结束、日志达到目标时间、没有NaN/Inf或致命错误，DAT完整。
02和04续算保持`firstprocess=.false.`；03自身中断后的续算需另行检查首次处理和时钟设置，
不能直接重复首次插绳命令。详细参数含义查看对应`.par.readme`，物理模型见[README.md](README.md)。

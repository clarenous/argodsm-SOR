package main

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"sync"
)

type SorTimes struct {
	Times    []float64 `json:"times"`
	max, min float64
	avg, mid float64
}

func (times *SorTimes) copy() []float64 {
	cpy := make([]float64, len(times.Times))
	copy(cpy, times.Times)
	return cpy
}

func (times *SorTimes) Max() float64 {
	if times.max != 0 {
		return times.max
	}

	if len(times.Times) == 0 {
		return 0
	}

	cpy := times.copy()
	sort.Float64s(cpy)
	times.max = cpy[len(cpy)-1]

	return times.max
}

func (times *SorTimes) Min() float64 {
	if times.min != 0 {
		return times.min
	}

	if len(times.Times) == 0 {
		return 0
	}

	cpy := times.copy()
	sort.Float64s(cpy)
	times.min = cpy[0]

	return times.min
}

func (times *SorTimes) Average() float64 {
	if times.avg != 0 {
		return times.avg
	}

	var sum float64
	for _, t := range times.Times {
		sum += t
	}
	times.avg = sum / float64(len(times.Times))

	return times.avg
}

func (times *SorTimes) Middle() float64 {
	if times.mid != 0 {
		return times.mid
	}

	if len(times.Times) == 0 {
		return 0
	}

	cpy := times.copy()
	sort.Float64s(cpy)
	if size := len(cpy); size%2 == 0 {
		times.mid = (cpy[size/2-1] + cpy[size/2]) / 2
	} else {
		times.mid = cpy[size/2]
	}

	return times.mid
}

type MethodType uint8

const (
	MethodRef MethodType = iota
	MethodMPI
	MethodArgoV1
	MethodArgoV2
	MethodArgoV3
	MethodEndTag
)

type Meta struct {
	Method    MethodType `json:"method"`
	Optimized bool       `json:"optimized"`
	NodeCount int        `json:"node_count"`
	DsmNX     int        `json:"dsm_nx"`
	DsmNY     int        `json:"dsm_ny"`
}

type ArraySize struct {
	I, J, K int
}

type StorageItem struct {
	Size  ArraySize `json:"size"`
	Meta  Meta      `json:"meta"`
	Times SorTimes  `json:"times"`
}

type StorageResult struct {
	Items []StorageItem `json:"items"`
}

type Result struct {
	mutex sync.RWMutex
	Items map[ArraySize]map[Meta]SorTimes `json:"items"`
}

func NewResult() *Result {
	return &Result{
		Items: map[ArraySize]map[Meta]SorTimes{},
	}
}

func NewResultFromBytes(data []byte) (*Result, error) {
	storage := &StorageResult{}
	if err := json.Unmarshal(data, storage); err != nil {
		return nil, err
	}
	result := NewResult()
	result.FromStorage(storage)
	return result, nil
}

func (result *Result) CSV() [][]string {
	storage := result.ToStorage()
	headers := []string{"im", "jm", "km", "method", "optimized", "node_count", "dsm_nx", "dsm_ny", "time"}
	var optimized = func(optimizedFlag bool) int {
		if optimizedFlag {
			return 1
		} else {
			return 0
		}
	}
	var itoa = func(i int) string {
		return strconv.Itoa(i)
	}
	var ftoa = func(f float64) string {
		return strconv.FormatFloat(f, 'f', 6, 64)
	}
	records := make([][]string, 0, 0)
	records = append(records, headers)
	for _, item := range storage.Items {
		record := []string{
			itoa(item.Size.I), itoa(item.Size.J), itoa(item.Size.K),
			itoa(int(item.Meta.Method)), itoa(optimized(item.Meta.Optimized)), itoa(item.Meta.NodeCount),
			itoa(item.Meta.DsmNX), itoa(item.Meta.DsmNY), ftoa(item.Times.Max()),
		}
		records = append(records, record)
	}
	return records
}

func (result *Result) Bytes() ([]byte, error) {
	return json.MarshalIndent(result.ToStorage(), "", "  ")
}

func (result *Result) ToStorage() *StorageResult {
	result.mutex.RLock()
	defer result.mutex.RUnlock()

	storage := &StorageResult{Items: []StorageItem{}}

	for size, metaMap := range result.Items {
		for meta, sorTimes := range metaMap {
			storage.Items = append(storage.Items, StorageItem{
				Size:  size,
				Meta:  meta,
				Times: sorTimes,
			})
		}
	}

	sort.SliceStable(storage.Items, func(i, j int) bool {
		x, y := storage.Items[i], storage.Items[j]

		if x.Size.I < y.Size.I {
			return true
		} else if x.Size.I > y.Size.I {
			return false
		}

		if x.Size.J < y.Size.J {
			return true
		} else if x.Size.J > y.Size.J {
			return false
		}

		if x.Size.K < y.Size.K {
			return true
		} else if x.Size.K > y.Size.K {
			return false
		}

		if x.Meta.Method < y.Meta.Method {
			return true
		} else if x.Meta.Method > y.Meta.Method {
			return false
		}

		if !x.Meta.Optimized && y.Meta.Optimized {
			return true
		} else if x.Meta.Optimized && !y.Meta.Optimized {
			return false
		}

		if x.Meta.NodeCount < y.Meta.NodeCount {
			return true
		} else if x.Meta.NodeCount > y.Meta.NodeCount {
			return false
		}

		if x.Meta.DsmNX < y.Meta.DsmNX {
			return true
		} else if x.Meta.DsmNX > y.Meta.DsmNX {
			return false
		}

		if x.Meta.DsmNY < y.Meta.DsmNY {
			return true
		} else if x.Meta.DsmNY > y.Meta.DsmNY {
			return false
		}

		return false
	})

	return storage
}

func (result *Result) FromStorage(storage *StorageResult) {
	result.mutex.Lock()
	defer result.mutex.Unlock()

	for i := range storage.Items {
		result.joinItem(storage.Items[i].Size, storage.Items[i].Meta, storage.Items[i].Times)
	}
}

func (result *Result) JoinItem(arraySize ArraySize, meta Meta, times SorTimes) {
	result.mutex.Lock()
	defer result.mutex.Unlock()

	result.joinItem(arraySize, meta, times)
}

func (result *Result) joinItem(arraySize ArraySize, meta Meta, times SorTimes) {
	se, ok := result.Items[arraySize]
	if !ok {
		se = map[Meta]SorTimes{}
		result.Items[arraySize] = se
	}
	se[meta] = times
}

func mustParseArraySize(sizeStr string) ArraySize {
	numStr := strings.Split(sizeStr, "_")
	i, _ := strconv.Atoi(numStr[0])
	j, _ := strconv.Atoi(numStr[1])
	k, _ := strconv.Atoi(numStr[2])
	return ArraySize{i, j, k}
}

func mustParseNumbers(filename string, numRegExp *regexp.Regexp) (numbers []int) {
	numStr := numRegExp.FindAllString(filename, -1)
	for _, str := range numStr {
		num, _ := strconv.Atoi(str)
		numbers = append(numbers, num)
	}
	return numbers
}

func collectFromRefFile(result *Result, path string, size ArraySize, optimize bool) {
	data, err := ioutil.ReadFile(path)
	if err != nil {
		log.Println("[ERROR]", size, optimize, path, err)
		return
	}

	const prefix = "SOR Time: "
	var times []float64

	lines := strings.Split(string(data), "\n")
	for _, line := range lines {
		if !strings.Contains(line, prefix) {
			continue
		}
		timeStr := strings.ReplaceAll(line, prefix, "")
		t, err := strconv.ParseFloat(timeStr, 64)
		if err != nil {
			log.Println("[ERROR]", size, optimize, timeStr, err)
			continue
		}
		times = append(times, t)
	}

	meta := Meta{Method: MethodRef, Optimized: optimize, NodeCount: 1}
	sorTimes := SorTimes{Times: times}

	result.JoinItem(size, meta, sorTimes)
}

func collectFromMpiFile(result *Result, path string, size ArraySize, optimize bool, metaNumbers []int) {
	data, err := ioutil.ReadFile(path)
	if err != nil {
		log.Println("[ERROR]", size, optimize, metaNumbers, path, err)
		return
	}

	const prefix = "SOR Time: "
	var times []float64

	lines := strings.Split(string(data), "\n")
	for _, line := range lines {
		if !strings.Contains(line, prefix) {
			continue
		}
		timeStr := strings.ReplaceAll(line, prefix, "")
		t, err := strconv.ParseFloat(timeStr, 64)
		if err != nil {
			log.Println("[ERROR]", size, optimize, metaNumbers, timeStr, err)
			continue
		}
		times = append(times, t)
	}

	meta := Meta{Method: MethodMPI, Optimized: optimize, NodeCount: metaNumbers[0]}
	sorTimes := SorTimes{Times: times}

	result.JoinItem(size, meta, sorTimes)
}

func collectFromArgoFile(result *Result, path string, size ArraySize, optimizationFlags [MethodArgoV3 - MethodArgoV1 + 1]bool, metaNumbers []int) {
	data, err := ioutil.ReadFile(path)
	if err != nil {
		log.Println("[ERROR]", size, optimizationFlags, metaNumbers, path, err)
		return
	}

	const prefix = "SOR: "
	var times []float64

	lines := strings.Split(string(data), "\n")
	for _, line := range lines {
		if !strings.Contains(line, prefix) {
			continue
		}
		timeStr := strings.ReplaceAll(line, prefix, "")
		t, err := strconv.ParseFloat(timeStr, 64)
		if err != nil {
			log.Println("[ERROR]", size, optimizationFlags, metaNumbers, timeStr, err)
			continue
		}
		times = append(times, t)
	}

	meta := Meta{
		Method:    MethodArgoV1 + MethodType(metaNumbers[0]-1),
		Optimized: optimizationFlags[metaNumbers[0]-1],
		NodeCount: metaNumbers[1],
		DsmNX:     metaNumbers[2],
		DsmNY:     metaNumbers[3],
	}
	sorTimes := SorTimes{Times: times}

	result.JoinItem(size, meta, sorTimes)
}

func collectFromBaseDir(result *Result, rootDir string, size ArraySize, optimizationFlags [MethodEndTag]bool) error {
	var refFileStr = "c.out"            // filename of c.out
	var mpiRegStr = `^mpi_[0-9]+\.out$` // filename like mpi_8.out
	mpiRegExp, err := regexp.Compile(mpiRegStr)
	if err != nil {
		return err
	}
	var argoRegStr = `^argo_[0-9]+_[0-9]+_[0-9]+_[0-9]+\.out$` // filename like argo_1_16_4_4.out
	argoRegExp, err := regexp.Compile(argoRegStr)
	if err != nil {
		return err
	}
	var numRegStr = `[0-9]+` // find numbers from filename
	numRegExp, err := regexp.Compile(numRegStr)
	if err != nil {
		return err
	}

	baseDir := filepath.Join(rootDir, fmt.Sprintf("%d_%d_%d", size.I, size.J, size.K))
	baseDirInfo, err := ioutil.ReadDir(baseDir)
	if err != nil {
		return err
	}
	for _, fileInfo := range baseDirInfo {
		if fileInfo.IsDir() {
			continue
		}
		if filename := fileInfo.Name(); filename == refFileStr {
			collectFromRefFile(result, filepath.Join(baseDir, filename), size, optimizationFlags[MethodRef])
		} else if mpiRegExp.MatchString(filename) {
			collectFromMpiFile(result, filepath.Join(baseDir, filename), size, optimizationFlags[MethodMPI],
				mustParseNumbers(filename, numRegExp))
		} else if argoRegExp.MatchString(filename) {
			collectFromArgoFile(result, filepath.Join(baseDir, filename), size,
				[MethodArgoV3 - MethodArgoV1 + 1]bool{optimizationFlags[MethodArgoV1], optimizationFlags[MethodArgoV2], optimizationFlags[MethodArgoV3]},
				mustParseNumbers(filename, numRegExp))
		}
	}

	return nil
}

func collectFromRootDir(result *Result, rootDir string, optimizationFlags [MethodEndTag]bool) error {
	var pathRegStr = `^[0-9]+_[0-9]+_[0-9]+$` // pathname like 512_512_256
	pathRegExp, err := regexp.Compile(pathRegStr)
	if err != nil {
		return err
	}

	// find all eligible paths
	var arraySizes []ArraySize
	targetDirInfo, err := ioutil.ReadDir(rootDir)
	if err != nil {
		return err
	}
	for _, pathInfo := range targetDirInfo {
		if pathRegExp.MatchString(pathInfo.Name()) {
			arraySizes = append(arraySizes, mustParseArraySize(pathInfo.Name()))
		}
	}

	// collect logs from each path
	for _, size := range arraySizes {
		if err := collectFromBaseDir(result, rootDir, size, optimizationFlags); err != nil {
			return err
		}
	}

	return nil
}

func main() {
	result := NewResult()
	if err := collectFromRootDir(result, "results", [MethodEndTag]bool{false, false, true, true, true}); err != nil {
		log.Fatalln("collectFromRootDir", "results", err)
	}
	if err := collectFromRootDir(result, "results-flags", [MethodEndTag]bool{true, true, false, false, false}); err != nil {
		log.Fatalln("collectFromRootDir", "results-flags", err)
	}

	data, err := result.Bytes()
	if err != nil {
		log.Fatalln("result.Bytes", err)
	}
	jsonFile, err := os.Create("test-results.json")
	if err != nil {
		log.Fatalln("os.Create", err)
	}
	defer jsonFile.Close()
	jsonFile.Write(data)

	csvFile, err := os.Create("test-results.csv")
	if err != nil {
		log.Fatalln("os.Create", err)
	}
	defer csvFile.Close()
	csvWriter := csv.NewWriter(csvFile)
	if err = csvWriter.WriteAll(result.CSV()); err != nil {
		log.Fatalln("csv.Writer.WriteAll", err)
	}
	csvWriter.Flush()
	csvFile.Sync()
}

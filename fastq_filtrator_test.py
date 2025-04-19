import os
import pytest
from fastq_filtrator import filter_fastq


@pytest.fixture
def create_input_fastq():
    file_path = "temp_input.fastq"
    content = ["@SRX079804\n", "ATGC\n", "+SRX079804\n", "DDDD\n"]
    with open(file_path, "w") as f:
        f.writelines(content)
    yield file_path
    try:
        os.remove(file_path)
    except PermissionError:
        pass


@pytest.fixture
def create_input_fastq_with_bed_quality():
    file_path = "temp_input_bed_quality.fastq"
    content = ["@SRX079804\n", "ATGC\n", "+SRX079804\n", "0000\n"]
    with open(file_path, "w") as f:
        f.writelines(content)
    yield file_path
    try:
        os.remove(file_path)
    except PermissionError:
        pass


@pytest.fixture
def create_input_fastq_with_low_gc():
    file_path = "temp_input_low_gc.fastq"
    content = ["@SRX079804\n", "AAAA\n", "+SRX079804\n", "DDDD\n"]
    with open(file_path, "w") as f:
        f.writelines(content)
    yield file_path
    try:
        os.remove(file_path)
    except PermissionError:
        pass


@pytest.fixture
def create_input_fastq_with_high_gc():
    file_path = "temp_input_high_gc.fastq"
    content = ["@SRX079804\n", "GCGC\n", "+SRX079804\n", "DDDD\n"]
    with open(file_path, "w") as f:
        f.writelines(content)
    yield file_path
    try:
        os.remove(file_path)
    except PermissionError:
        pass


@pytest.fixture
def create_input_fastq_with_low_length():
    file_path = "temp_input_low_length.fastq"
    content = ["@SRX079804\n", "AGCT\n", "+SRX079804\n", "DDDD\n"]
    with open(file_path, "w") as f:
        f.writelines(content)
    yield file_path
    try:
        os.remove(file_path)
    except PermissionError:
        pass


@pytest.fixture
def create_input_fastq_with_high_length():
    file_path = "temp_input_high_length.fastq"
    content = [
        "@SRX079804\n",
        "AGCTAGCTAGCTAGCT\n",
        "+SRX079804\n",
        "DDDDDDDDDDDDDDDD\n",
    ]
    with open(file_path, "w") as f:
        f.writelines(content)
    yield file_path
    try:
        os.remove(file_path)
    except PermissionError:
        pass


@pytest.fixture
def tmp_file():
    file_path = "tmp_filtered.fastq"
    if os.path.exists(file_path):
        os.remove(file_path)
    yield file_path
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
        except PermissionError:
            pass


class TestWriteFasta:
    def test_filter_fastq_exists(self, create_input_fastq, tmp_file):
        filter_fastq(input_fastq=create_input_fastq, output_fastq=tmp_file)
        assert os.path.exists(tmp_file)

    def test_filter_fastq_valid_content(self, create_input_fastq, tmp_file):
        filter_fastq(input_fastq=create_input_fastq, output_fastq=tmp_file)
        with open(tmp_file, "r") as f:
            lines = f.readlines()
        assert len(lines) % 4 == 0


class TestFastqFilter:
    def test_low_quality_fails(self, create_input_fastq_with_bed_quality, tmp_file):
        filter_fastq(
            input_fastq=create_input_fastq_with_bed_quality,
            quality_threshold=20,
            output_fastq=tmp_file,
        )

        with open(tmp_file, "r") as f:
            lines = f.readlines()

        assert len(lines) == 0

    def test_gc_fails_low(self, create_input_fastq_with_low_gc, tmp_file):
        filter_fastq(
            input_fastq=create_input_fastq_with_low_gc,
            gc_bounds=(40, 60),
            output_fastq=tmp_file,
        )

        with open(tmp_file, "r") as f:
            lines = f.readlines()

        assert len(lines) == 0

    def test_gc_fails_high(self, create_input_fastq_with_high_gc, tmp_file):
        filter_fastq(
            input_fastq=create_input_fastq_with_high_gc,
            gc_bounds=(40, 60),
            output_fastq=tmp_file,
        )

        with open(tmp_file, "r") as f:
            lines = f.readlines()

        assert len(lines) == 0

    def test_length_fails_low(self, create_input_fastq_with_low_length, tmp_file):
        filter_fastq(
            input_fastq=create_input_fastq_with_low_length,
            length_bounds=(5, 100),
            output_fastq=tmp_file,
        )

        with open(tmp_file, "r") as f:
            lines = f.readlines()

        assert len(lines) == 0

    def test_length_fails_high(self, create_input_fastq_with_high_length, tmp_file):
        filter_fastq(
            input_fastq=create_input_fastq_with_high_length,
            length_bounds=(0, 15),
            output_fastq=tmp_file,
        )

        with open(tmp_file, "r") as f:
            lines = f.readlines()

        assert len(lines) == 0


class TestError:
    def test_error(self, create_input_fastq, tmp_file):
        with pytest.raises(ValueError):
            filter_fastq(
                input_fastq=create_input_fastq,
                length_bounds=(100, 0),
                output_fastq=tmp_file,
            )

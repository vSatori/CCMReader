#include "vtkCCMReader.h"
#include <map>
#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkDataSetAttributes.h>
#include "ccmio.h"
#include "ccmiotypes.h"


vtkStandardNewMacro(vtkCCMReader)

static void read_map(CCMIOID map, std::vector<int>& map_data)
{
    CCMIOSize size = 0;
    CCMIOEntitySize(nullptr, map, &size, nullptr);
    map_data.resize(size);
    CCMIOReadMap(nullptr, map, &map_data[0], kCCMIOStart, kCCMIOEnd);
    return;
}


static void read_faces(CCMIOID face, CCMIOEntity type, std::map<int, std::vector<int>>& cell_faces, std::map<int, std::vector<int>>& face_nodes)
{
    CCMIOSize size = 0;
    std::vector<int> face_cells;
    CCMIOEntitySize(nullptr, face, &size, nullptr);
    face_cells.resize(type == CCMIOEntity::kCCMIOInternalFaces ? size * 2 : size);
    CCMIOReadFaceCells(nullptr, face, type, &face_cells[0], kCCMIOStart, kCCMIOEnd);
    CCMIOID map;
    CCMIOSize streamSize = 0;
    std::vector<int> stream_data;
    CCMIOReadFaces(nullptr, face, type, &map, &streamSize, nullptr, kCCMIOStart, kCCMIOEnd);
    stream_data.resize(streamSize);
    CCMIOReadFaces(nullptr, face, type, &map, &streamSize, &stream_data[0], kCCMIOStart, kCCMIOEnd);
    std::vector<int> map_data;
    read_map(map, map_data);

    if (type == CCMIOEntity::kCCMIOBoundaryFaces)
    {
        for (int index = 0; index < face_cells.size(); ++index)
        {
            auto cell = face_cells[index];
            if (cell_faces.find(cell) != cell_faces.end())
            {
                auto& cell_face = cell_faces[cell];
                cell_face.push_back(map_data[index]);
            }
            else
            {
                std::vector<int> cell_face;
                cell_face.push_back(map_data[index]);
                cell_faces[cell] = std::move(cell_face);
            }
        }
    }
    else
    {
        std::map<int, std::vector<int>> nbr_ids;
        for (int index = 0; index < face_cells.size(); ++index)
        {
            if (index % 2)
            {
                continue;
            }
            auto cell = face_cells[index];
            if (nbr_ids.find(cell) != nbr_ids.end())
            {
                auto& faces = nbr_ids[cell];
                faces.push_back(map_data[index / 2]);
            }
            else
            {
                std::vector<int> faces;
                faces.push_back(map_data[index / 2]);
                nbr_ids[cell] = std::move(faces);
            }

        }

        for (int index = 0; index < face_cells.size(); ++index)
        {
            if (index % 2 == 0)
            {
                continue;
            }
            auto cell = face_cells[index];
            if (cell_faces.find(cell) != cell_faces.end())
            {
                auto& cell_face = cell_faces[cell];
                cell_face.push_back(map_data[index / 2]);

            }
            else
            {
                std::vector<int> cell_face;
                cell_face.push_back(map_data[index / 2]);
                if (nbr_ids.find(cell) != nbr_ids.end())
                {
                    auto& faces = nbr_ids[cell];
                    for (const auto& face : faces)
                    {
                        cell_face.push_back(face);
                    }
                }
                cell_faces[cell] = std::move(cell_face);
            }
        }

        for (auto it = nbr_ids.begin(); it != nbr_ids.end(); ++it)
        {
            if (cell_faces.find(it->first) == cell_faces.end())
            {
                cell_faces[it->first] = std::move(it->second);
            }
        }
    }

    int current_index = 0;
    int index_count = stream_data[0];
    current_index += 1;
    int face_index = 0;
    while (true)
    {
        std::vector<int> indices;
        indices.reserve(index_count);
        int end = current_index + index_count;
        for (; current_index < end; ++current_index)
        {
            indices.push_back(stream_data[current_index]);
        }
        face_nodes[map_data[face_index]] = indices;
        face_index += 1;
        if (current_index >= stream_data.size() - 1)
        {
            break;
        }
        index_count = stream_data[current_index];
        current_index += 1;
    }

}



static void read_field(CCMIOID state, vtkUnstructuredGrid* grid, const std::map<int, int>& vertex_index_map, const std::map<int, int>& cell_index_map)
{
    CCMIOID field_set;
    int i = 0;
    while (CCMIOError::kCCMIONoErr == CCMIONextEntity(nullptr, state, kCCMIOFieldSet, &i, &field_set))
    {
        CCMIOID field_phase;
        CCMIOID field;
        int j = 0;
        if (CCMIOError::kCCMIONoErr != CCMIONextEntity(nullptr, field_set, kCCMIOFieldPhase, &j, &field_phase))
        {
            break;
        }
        j = 0;
        while (CCMIOError::kCCMIONoErr == CCMIONextEntity(nullptr, field_phase, kCCMIOField, &j, &field))
        {
            char name[64];
            char short_name[64];
            CCMIOSize size = 0;
            CCMIODimensionality field_dim;
            CCMIODataType dt;
            CCMIOID map;
            CCMIODataLocation loc;
            CCMIOID field_data_id;
            int count = 0;
            CCMIOReadField(nullptr, field, name, short_name, &field_dim, &dt);
            if (dt > 3)
            {
                continue;
            }
            CCMIODataLocation current_loc = CCMIODataLocation::kCCMIOCell;
            vtkDataSetAttributes* attr = nullptr;
            vtkNew<vtkDoubleArray> arr;
            arr->SetName(name);
            arr->SetNumberOfComponents(field_dim);

            int k = 0;
            std::map<int, double> field_data_map;
            while (CCMIOError::kCCMIONoErr == CCMIONextEntity(nullptr, field, CCMIOEntity::kCCMIOFieldData, &k, &field_data_id))
            {
                CCMIOEntitySize(nullptr, field_data_id, &size, nullptr);
                int factor = 1;
                if (field_dim == CCMIODimensionality::kCCMIOVector) { factor = 3; }
                if (field_dim == CCMIODimensionality::kCCMIOTensor) { factor = 9; }

                if (dt == CCMIODataType::kCCMIOInt32 || dt == CCMIODataType::kCCMIOInt64)
                {
                    break;
                    std::vector<int> field_data;
                    field_data.resize(size * factor);
                    CCMIOReadFieldDatai(nullptr, field_data_id, &map, &loc, &field_data[0], kCCMIOStart, kCCMIOEnd);
                }
                else if (dt == CCMIODataType::kCCMIOFloat32)
                {
                    std::vector<float> field_data;
                    field_data.resize(size * factor);
                    CCMIOReadFieldDataf(nullptr, field_data_id, &map, &loc, &field_data[0], kCCMIOStart, kCCMIOEnd);
                    if (loc > 1)
                    {
                        continue;
                    }
                    if (attr == nullptr)
                    {
                        attr = (loc == CCMIODataLocation::kCCMIOVertex) ? grid->GetAttributes(vtkDataSet::POINT) : grid->GetAttributes(vtkDataSet::CELL);
                    }
                    std::vector<int> field_map_data;
                    read_map(map, field_map_data);
                    for (int index = 0; index < size; ++index)
                    {
                        field_data_map[field_map_data[index]] = field_data[index];
                    }
                }
                else if (dt == CCMIODataType::kCCMIOFloat64)
                {
                    std::vector<double> field_data;
                    field_data.resize(size * factor);
                    CCMIOReadFieldDatad(nullptr, field_data_id, &map, &loc, &field_data[0], kCCMIOStart, kCCMIOEnd);
                    if (loc > 1)
                    {
                        continue;
                    }
                    if (attr == nullptr)
                    {
                        current_loc = loc;
                        attr = (loc == CCMIODataLocation::kCCMIOVertex) ? grid->GetAttributes(vtkDataSet::POINT) : grid->GetAttributes(vtkDataSet::CELL);
                    }
                    std::vector<int> field_map_data;
                    read_map(map, field_map_data);
                    for (int index = 0; index < size; ++index)
                    {
                        field_data_map[field_map_data[index]] = field_data[index];
                    }
                }
            }
            if (attr)
            {
                const std::map<int, int>* id_map = nullptr;
                if (current_loc == CCMIODataLocation::kCCMIOCell)
                {
                    id_map = &cell_index_map;
                    arr->SetNumberOfTuples(grid->GetNumberOfCells());
                }
                else
                {
                    id_map = &vertex_index_map;
                    arr->SetNumberOfTuples(grid->GetNumberOfPoints());
                }
                arr->Fill(0.0);
                for (auto it = field_data_map.begin(); it != field_data_map.end(); ++it)
                {
                    int cell_index = it->first;
                    if (cell_index_map.find(cell_index) == cell_index_map.end())
                    {
                        continue;
                    }
                    vtkIdType id = cell_index_map.at(cell_index);
                    arr->SetTuple1(id, it->second);
                }
                attr->AddArray(arr);
            }


        }
    }
}

static bool read_grid(CCMIOID state_id, vtkUnstructuredGrid* grid, std::map<int, int>& vertex_index_map, std::map<int, int>& cell_index_map)
{
    CCMIOID vertex_id;
    CCMIOID topo_id;
    CCMIOID processor_id;
    int id_index = 0;
    int dims = 0;
    float scale = 1.f;
    std::vector<float> vertices_data;
    CCMIOSize points_count = 0;

    CCMIOID internal_face;
    CCMIOID boundary_face;
    CCMIOID cells;
    CCMIOID cell_map;
    CCMIOSize cells_count = 0;
    CCMIOSize stream_size = 0;
    CCMIOSize points_map_data_count = 0;
    CCMIOID points_map;
    std::vector<int> points_map_data;

    id_index = 0;
    CCMIONextEntity(nullptr, state_id, CCMIOEntity::kCCMIOProcessor, &id_index, &processor_id);
    CCMIOReadProcessor(nullptr, processor_id, &vertex_id, &topo_id, nullptr, nullptr);
    if (topo_id.type == CCMIOEntity::kCCMIOMaxEntity || vertex_id.type == CCMIOEntity::kCCMIOMaxEntity)
    {
        return false;
    }

    //read vertices
    id_index = 0;
    CCMIONextEntity(nullptr, topo_id, CCMIOEntity::kCCMIOVertices, &id_index, &vertex_id);
    CCMIOEntitySize(nullptr, vertex_id, &points_count, nullptr);
    CCMIOReadVerticesf(nullptr, vertex_id, &dims, &scale, &points_map, nullptr, kCCMIOStart, kCCMIOEnd);
    vertices_data.resize(dims * points_count);
    CCMIOReadVerticesf(nullptr, vertex_id, &dims, &scale, &points_map, &vertices_data[0], kCCMIOStart, kCCMIOEnd);
    read_map(points_map, points_map_data);

    for (int index = 0; index < points_map_data.size(); ++index)
    {
        vertex_index_map[points_map_data[index]] = index;
    }

    //read cells
    id_index = 0;
    CCMIONextEntity(nullptr, topo_id, CCMIOEntity::kCCMIOCells, &id_index, &cells);
    CCMIOEntitySize(nullptr, cells, &cells_count, nullptr);
    std::vector<int> cell_types;
    cell_types.resize(cells_count);
    CCMIOReadCells(nullptr, cells, &cell_map, &cell_types[0], kCCMIOStart, kCCMIOEnd);
    std::vector<int> cell_map_data;
    read_map(cell_map, cell_map_data);
    //read faces
    std::map<int, std::vector<int>> cell_faces;
    std::map<int, std::vector<int>> face_nodes;
    id_index = 0;
    CCMIONextEntity(nullptr, topo_id, CCMIOEntity::kCCMIOInternalFaces, &id_index, &internal_face);
    read_faces(internal_face, CCMIOEntity::kCCMIOInternalFaces, cell_faces, face_nodes);
    id_index = 0;
    while (CCMIONextEntity(nullptr, topo_id, CCMIOEntity::kCCMIOBoundaryFaces, &id_index, &boundary_face) == CCMIOError::kCCMIONoErr)
    {
        read_faces(boundary_face, CCMIOEntity::kCCMIOBoundaryFaces, cell_faces, face_nodes);
    }

    //make grid
    vtkNew<vtkPoints> points;
    points->Allocate(100000);
    for (int pIndex = 0; pIndex < points_count; ++pIndex)
    {
        auto x = vertices_data[pIndex * dims + 0];
        auto y = vertices_data[pIndex * dims + 1];
        double z = 0.0;
        if (dims == 3)
        {
            z = vertices_data[pIndex * 3 + 2];
        }
        points->InsertNextPoint(x, y, z);

    }
    grid->SetPoints(points);
    grid->Allocate();

    for (auto it = cell_faces.begin(); it != cell_faces.end(); ++it)
    {
        const auto& cell_face = it->second;
        std::vector<vtkIdType> ids;
        for (int index = 0; index < cell_face.size(); ++index)
        {
            const auto& face_index = cell_face[index];
            const auto& face = face_nodes[face_index];
            ids.push_back(face.size());
            for (const auto& node_index : face)
            {
                ids.push_back(vertex_index_map[node_index]);
            }
        }
        cell_index_map[it->first] = grid->InsertNextCell(VTKCellType::VTK_POLYHEDRON, cell_face.size(), &ids[0]);
    }
    return true;
}


void vtkCCMReader::PrintSelf(ostream & os, vtkIndent indent)
{
    Superclass::PrintSelf(os, indent);
}

vtkCCMReader::vtkCCMReader()
{
    this->FileName = nullptr;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

vtkCCMReader::~vtkCCMReader()
{
    if (this->FileName)
    {
        delete this->FileName;
    }
}

int vtkCCMReader::RequestData(vtkInformation * request, vtkInformationVector ** inputVector, vtkInformationVector * outputVector)
{
    auto outInfo = outputVector->GetInformationObject(0);
    auto grid = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    if (!this->FileName)
    {
        return VTK_ERROR;
    }
    CCMIOID root_id;
    if (CCMIOOpenFile(nullptr, this->FileName, CCMIOIOType::kCCMIORead, &root_id) != CCMIOError::kCCMIONoErr)
    {
        return VTK_ERROR;
    }
    CCMIOID state_id;
    int id_index = 0;
    CCMIONextEntity(nullptr, root_id, CCMIOEntity::kCCMIOState, &id_index, &state_id);
    std::map<int, int> vertex_index_map;
    std::map<int, int> cell_index_map;
    if (!read_grid(state_id, grid, vertex_index_map, cell_index_map))
    {
        CCMIOCloseFile(nullptr, root_id);
        return VTK_ERROR;
    }
    read_field(state_id, grid, vertex_index_map, cell_index_map);
    CCMIOCloseFile(nullptr, root_id);
    return VTK_OK;
}

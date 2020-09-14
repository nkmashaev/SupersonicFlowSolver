#include <memory>
#include <vector>
#include <map>
#ifndef _PSST_BoundaryCondition_H_
#define _PSST_BoundaryCondition_H_

struct psstVolumeStorage;
class psstInterfaceBoundaryCondition
{
private:

public:
	virtual void initialize() = 0;
	virtual psstInterfaceBoundaryCondition * clone() const = 0;

	virtual void calculate_at(size_t i) const = 0;
	virtual void get_data(psstVolumeStorage * storage) = 0;
	virtual size_t type() const = 0;
	virtual ~psstInterfaceBoundaryCondition() {};
};

struct psstBoundary
{
public:
	std::shared_ptr<psstInterfaceBoundaryCondition> dynamic_cond;
	std::shared_ptr<psstInterfaceBoundaryCondition> thermal_cond;
};

class psstSupersonicInlet : public psstInterfaceBoundaryCondition
{
private:
	std::vector<double> velocity_;
	double pressure_;
	psstVolumeStorage * storage_ = 0;

	enum dyn_bc_id
	{
		ss_inlet = 23,	//������������� ������� �������
		ss_outlet = 26,	//������������� �������� �������
		wall_no_slip = 15,	//������ � ����������������
		inlet = 27,	//���������� �������� �������
		outlet = 25	//���������� �������� �������
	};

	psstSupersonicInlet() {};
public:
	psstSupersonicInlet(const std::vector<double> & velocity, double pressure);

	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const {};
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return ss_inlet; };
	virtual ~psstSupersonicInlet() {};
};

class psstSupersonicOutlet : public psstInterfaceBoundaryCondition
{
private:
	psstVolumeStorage * storage_ = 0;

	enum dyn_bc_id
	{
		ss_inlet = 23,	//������������� ������� �������
		ss_outlet = 26,	//������������� �������� �������
		wall_no_slip = 15,	//������ � ����������������
		inlet = 27,	//���������� �������� �������
		outlet = 25	//���������� �������� �������
	};

public:
	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const;
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return ss_outlet; };
	virtual ~psstSupersonicOutlet() {};
};

class psstInviscidWall : public psstInterfaceBoundaryCondition
{
private:
	psstVolumeStorage * storage_ = 0;

	enum dyn_bc_id
	{
		ss_inlet = 23,	//������������� ������� �������
		ss_outlet = 26,	//������������� �������� �������
		wall_no_slip = 15,	//������ � ����������������
		inlet = 27,	//���������� �������� �������
		outlet = 25	//���������� �������� �������
	};

public:
	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const;
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return wall_no_slip; };
	virtual ~psstInviscidWall() {};
};

class psstInlet : public psstInterfaceBoundaryCondition
{
private:
	psstVolumeStorage * storage_ = 0;
	std::vector<double> velocity_;

	enum dyn_bc_id
	{
		ss_inlet = 23,	//������������� ������� �������
		ss_outlet = 26,	//������������� �������� �������
		wall_no_slip = 15,	//������ � ����������������
		inlet = 27,	//���������� �������� �������
		outlet = 25	//���������� �������� �������
	};

	psstInlet();
public:
	psstInlet(const std::vector<double> &);
	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const;
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return inlet; };
	virtual ~psstInlet() {};
};

class psstOutlet : public psstInterfaceBoundaryCondition
{
private:
	psstVolumeStorage * storage_ = 0;
	double pressure_;

	enum dyn_bc_id
	{
		ss_inlet = 23,	//������������� ������� �������
		ss_outlet = 26,	//������������� �������� �������
		wall_no_slip = 15,	//������ � ����������������
		inlet = 27,	//���������� �������� �������
		outlet = 25	//���������� �������� �������
	};

public:
	psstOutlet(double pressure);
	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const;
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return outlet; };
	virtual ~psstOutlet() {};
};

//"��������" �������

class psstThermalInlet : public psstInterfaceBoundaryCondition
{
private:
	psstVolumeStorage * storage_ = 0;
	double temperature_;

	enum term_bcID
	{
		wall_T = 11,	//������. T ��������
		wall_Q = 18,	//������. �������� �������� �����
		thermal_inlet = 27,	//������� �������
		thermal_outlet = 25	//�������� �������
	};

	psstThermalInlet();
public:
	psstThermalInlet(double temparature);
	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const {};
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return thermal_inlet; };
	virtual ~psstThermalInlet() {};
};

class psstThermalOutlet : public psstInterfaceBoundaryCondition
{
private:
	psstVolumeStorage * storage_ = 0;

	enum term_bcID
	{
		wall_T = 11,	//������. T ��������
		wall_Q = 18,	//������. �������� �������� �����
		thermal_inlet = 27,	//������� �������
		thermal_outlet = 25	//�������� �������
	};

public:
	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const;
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return thermal_outlet; };
	virtual ~psstThermalOutlet() {};
};

class psstThermalWall_T : public psstInterfaceBoundaryCondition
{
private:
	psstVolumeStorage * storage_ = 0;
	double temperature_;

	enum term_bcID
	{
		wall_T = 11,	//������. T ��������
		wall_Q = 18,	//������. �������� �������� �����
		thermal_inlet = 27,	//������� �������
		thermal_outlet = 25	//�������� �������
	};

	psstThermalWall_T();
public:

	psstThermalWall_T(double temperature);
	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const {};
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return wall_T; };
	virtual ~psstThermalWall_T() {};
};

class psstThermalWall_Q : public psstInterfaceBoundaryCondition
{
private:
	psstVolumeStorage * storage_ = 0;
	enum term_bcID
	{
		wall_T = 11,	//������. T ��������
		wall_Q = 18,	//������. �������� �������� �����
		thermal_inlet = 27,	//������� �������
		thermal_outlet = 25	//�������� �������
	};

public:
	virtual void initialize();
	virtual psstInterfaceBoundaryCondition * clone() const;

	virtual void calculate_at(size_t i) const;
	virtual void get_data(psstVolumeStorage * storage);
	virtual size_t type() const { return wall_Q; };
	virtual ~psstThermalWall_Q() {};
};
#endif
